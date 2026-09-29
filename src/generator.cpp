#include "pvmls/generator.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

extern "C" {
#include <jpeglib.h>
}

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace pvmls {
namespace {

namespace fs = std::filesystem;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Itag = CGAL::Exact_predicates_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>;
using Point = K::Point_2;
using VH = CDT::Vertex_handle;

constexpr double kGeomTol = 2e-11;
constexpr double kAreaTol = 2e-8;
constexpr const char* kSchemaVersion = "1.0.0";

struct P2 { double x{}, y{}; };
struct InterfaceGeom {
    P2 a, b;
    double pa{-1.0}, pb{-1.0};
    int edge_a{-1}, edge_b{-1};
    double sign_factor{1.0};
};
struct Node {
    int id{-1};
    P2 p;
    double phi{};
    int sign{};
    std::string constraint{"interior"};
    int parent_edge{-1};
    int parent_segment{-1};
    double parameter{-1.0};
    bool primary{true};
    bool fixed{false};
};
struct Tri {
    std::array<int,3> v{};
    int phase{};
    double quality{};
    std::array<int,6> p2{};
    std::array<int,3> pressure{};
};
struct Quad {
    std::array<int,4> v{};
    int phase{};
    double quality{};
    std::array<int,9> q2{};
    std::array<int,4> pressure{};
};
struct PressureRecord { int id{}, geom{}, phase{}; };
struct Mesh {
    std::vector<Node> nodes;
    std::vector<Tri> tris;
    std::vector<Quad> quads;
    std::vector<PressureRecord> pressures;
    std::vector<std::array<int,2>> boundary_edges;
    std::vector<std::array<int,2>> interface_edges;
    InterfaceGeom iface;
    std::array<double,4> corner_phi{};
    std::array<int,4> corner_sign{};
    double tri_min{1.0}, tri_mean{1.0}, quad_min{1.0}, quad_mean{1.0};
    double quad_count_fraction{}, quad_area_fraction{};
    std::array<double,2> phase_area{};
};
struct Candidate {
    Mesh mesh;
    std::string signature;
};

const std::array<P2,4> corners{{{-1,-1},{1,-1},{1,1},{-1,1}}};

double cross(P2 a, P2 b, P2 c) {
    return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
}
double dist2(P2 a, P2 b) {
    const double dx=a.x-b.x, dy=a.y-b.y;
    return dx*dx+dy*dy;
}
double area_poly(const std::vector<P2>& p) {
    double a=0.0;
    for (std::size_t i=0;i<p.size();++i) {
        const auto& u=p[i]; const auto& v=p[(i+1)%p.size()];
        a += u.x*v.y-u.y*v.x;
    }
    return 0.5*a;
}
double phi_raw(const InterfaceGeom& g, P2 p) {
    const double dx=g.b.x-g.a.x, dy=g.b.y-g.a.y;
    const double n=std::hypot(dx,dy);
    if (!(n>0.0)) throw std::runtime_error("degenerate interface");
    return g.sign_factor*((dx*(p.y-g.a.y)-dy*(p.x-g.a.x))/n);
}
std::array<int,4> expected_signs(int cid) {
    if (cid==2) return {-1,1,1,1};
    if (cid==3) return {-1,-1,1,1};
    if (cid==6) return {0,-1,1,1};
    if (cid==11) return {0,1,0,-1};
    throw std::invalid_argument("unsupported case");
}
P2 edge_point(int edge, double t) {
    const auto a=corners[static_cast<std::size_t>(edge)];
    const auto b=corners[static_cast<std::size_t>((edge+1)%4)];
    return {a.x+t*(b.x-a.x), a.y+t*(b.y-a.y)};
}
double edge_parameter(int edge, P2 p) {
    const auto a=corners[static_cast<std::size_t>(edge)];
    const auto b=corners[static_cast<std::size_t>((edge+1)%4)];
    const double dx=b.x-a.x, dy=b.y-a.y;
    return ((p.x-a.x)*dx+(p.y-a.y)*dy)/(dx*dx+dy*dy);
}
int boundary_edge(P2 p) {
    if (std::abs(p.y+1.0)<=kGeomTol && p.x>=-1-kGeomTol && p.x<=1+kGeomTol) return 0;
    if (std::abs(p.x-1.0)<=kGeomTol && p.y>=-1-kGeomTol && p.y<=1+kGeomTol) return 1;
    if (std::abs(p.y-1.0)<=kGeomTol && p.x>=-1-kGeomTol && p.x<=1+kGeomTol) return 2;
    if (std::abs(p.x+1.0)<=kGeomTol && p.y>=-1-kGeomTol && p.y<=1+kGeomTol) return 3;
    return -1;
}
bool at_corner(P2 p, int* out=nullptr) {
    for (int i=0;i<4;++i) if (dist2(p,corners[static_cast<std::size_t>(i)])<1e-24) {
        if (out) *out=i; return true;
    }
    return false;
}
double triangle_quality(P2 a,P2 b,P2 c) {
    const double A=0.5*cross(a,b,c);
    if (!(A>0.0)) return -1.0;
    const double s=dist2(a,b)+dist2(b,c)+dist2(c,a);
    return 4.0*std::sqrt(3.0)*A/s;
}
double quad_quality(const std::array<P2,4>& p) {
    double q=std::numeric_limits<double>::infinity();
    for (int i=0;i<4;++i) {
        const auto& c=p[static_cast<std::size_t>(i)];
        const auto& n=p[static_cast<std::size_t>((i+1)%4)];
        const auto& prev=p[static_cast<std::size_t>((i+3)%4)];
        const double ax=n.x-c.x, ay=n.y-c.y;
        const double bx=prev.x-c.x, by=prev.y-c.y;
        const double den=std::hypot(ax,ay)*std::hypot(bx,by);
        if (!(den>0.0)) return -1.0;
        const double sj=(ax*by-ay*bx)/den;
        q=std::min(q,sj);
    }
    return q;
}
double quad_area(const std::array<P2,4>& p) {
    return std::abs(area_poly(std::vector<P2>(p.begin(),p.end())));
}
std::uint64_t mix_seed(std::uint64_t s, int cid, int tid, int attempt) {
    auto mix=[](std::uint64_t x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x^(x>>30))*0xbf58476d1ce4e5b9ULL;
        x = (x^(x>>27))*0x94d049bb133111ebULL;
        return x^(x>>31);
    };
    return mix(s ^ mix(static_cast<std::uint64_t>(cid)*1315423911ULL)
                 ^ mix(static_cast<std::uint64_t>(tid)*2654435761ULL)
                 ^ mix(static_cast<std::uint64_t>(attempt+1)));
}
double stratified(int idx, int n, double lo, double hi, std::mt19937_64& rng, int salt) {
    if (n<=1) return 0.5*(lo+hi);
    if (idx==0) return lo;
    if (idx==n-1) return hi;
    std::uniform_real_distribution<double> U(0.15,0.85);
    const double u=(static_cast<double>(idx)+std::fmod(U(rng)+0.173*salt,1.0))/static_cast<double>(n);
    return lo+(hi-lo)*std::clamp(u,0.0,1.0);
}
InterfaceGeom sample_interface(const Config& cfg,int cid,int tid,int attempt,std::mt19937_64& rng) {
    const double lo=cfg.min_edge_fraction, hi=1.0-cfg.min_edge_fraction;
    InterfaceGeom g;
    if (cid==2) {
        const double t0=stratified(tid,cfg.templates_per_case,lo,hi,rng,1);
        const double t3=1.0-stratified(tid,cfg.templates_per_case,lo,hi,rng,2);
        g={edge_point(0,t0),edge_point(3,t3),t0,t3,0,3,1.0};
    } else if (cid==3) {
        const double t1=stratified(tid,cfg.templates_per_case,lo,hi,rng,3);
        const int shifted=(tid+std::max(1,cfg.templates_per_case/2))%cfg.templates_per_case;
        const double t3=stratified(shifted,cfg.templates_per_case,lo,hi,rng,4);
        g={edge_point(1,t1),edge_point(3,t3),t1,t3,1,3,1.0};
    } else if (cid==6) {
        double t1=stratified(tid,cfg.templates_per_case,lo,hi,rng,5);
        if(attempt>0){
            // At corner 0 the negative phase contains a wedge with angle atan(t1).
            // Its best possible conforming triangle quality is bounded by that angle.
            // Preserve the exact near-limit sample as attempt 0, but on rejection
            // deterministically resample above the quality-admissible lower bound.
            auto wedge_quality=[](double t){
                const double theta=std::atan(t);
                return std::sqrt(3.0)*std::sin(theta)/(2.0-std::cos(theta));
            };
            double a=lo,b=hi;
            if(wedge_quality(b)<cfg.min_triangle_quality)
                throw std::runtime_error("case 6 quality threshold is incompatible with permitted interface interval");
            for(int k=0;k<60;++k){
                const double mid=0.5*(a+b);
                if(wedge_quality(mid)<cfg.min_triangle_quality)a=mid;else b=mid;
            }
            const double safe=std::min(hi,b*1.15+1e-6);
            t1=std::max(t1,safe);
        }
        g={corners[0],edge_point(1,t1),0.0,t1,0,1,1.0};
    } else if (cid==11) {
        g={corners[0],corners[2],0.0,1.0,0,2,1.0};
    } else throw std::invalid_argument("unsupported case");

    const auto want=expected_signs(cid);
    auto score=[&](double sf){
        g.sign_factor=sf; int ok=0;
        for (int i=0;i<4;++i) {
            const int s=classify_phi(phi_raw(g,corners[static_cast<std::size_t>(i)]),cfg.phi_zero_tol);
            if (s==want[static_cast<std::size_t>(i)]) ++ok;
        }
        return ok;
    };
    if (score(1.0)!=4) {
        if (score(-1.0)!=4) throw std::runtime_error("sampled interface does not realize requested sign family");
    }
    return g;
}
Node make_node(P2 p,const InterfaceGeom& g,const Config& cfg,bool primary=true) {
    Node n; n.p=p; n.primary=primary;
    int ci=-1;
    const int be=boundary_edge(p);
    const double pr=phi_raw(g,p);
    const bool on_iface=std::abs(pr)<=5e-10;
    if (at_corner(p,&ci)) {
        n.constraint="corner"; n.parent_edge=ci; n.fixed=true;
    } else if (on_iface) {
        n.constraint="interface"; n.parent_segment=0;
        n.fixed=(dist2(p,g.a)<=1e-22 || dist2(p,g.b)<=1e-22);
        n.parameter=((p.x-g.a.x)*(g.b.x-g.a.x)+(p.y-g.a.y)*(g.b.y-g.a.y))/dist2(g.a,g.b);
        n.parent_edge=be;
    } else if (be>=0) {
        n.constraint="boundary"; n.parent_edge=be; n.parameter=edge_parameter(be,p);
    } else {
        n.constraint="interior";
    }
    if (on_iface) n.phi=0.0;
    else n.phi=pr;
    n.sign=classify_phi(n.phi,cfg.phi_zero_tol);
    return n;
}
std::array<int,2> edge_key(int a,int b) { return a<b?std::array<int,2>{a,b}:std::array<int,2>{b,a}; }

void add_constraint_chain(CDT& cdt,const std::vector<VH>& chain) {
    for (std::size_t i=1;i<chain.size();++i) cdt.insert_constraint(chain[i-1],chain[i]);
}
std::vector<P2> sample_segment(P2 a,P2 b,double h,bool include_ends=true) {
    const int n=std::max(1,static_cast<int>(std::ceil(std::sqrt(dist2(a,b))/h)));
    std::vector<P2> out;
    const int s=include_ends?0:1, e=include_ends?n:n-1;
    for (int i=s;i<=e;++i) {
        const double t=static_cast<double>(i)/n;
        out.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});
    }
    return out;
}
bool inside_square(P2 p) {
    return p.x>=-1.0-kGeomTol&&p.x<=1.0+kGeomTol&&p.y>=-1.0-kGeomTol&&p.y<=1.0+kGeomTol;
}

Mesh triangulate_primary(const Config& cfg,int cid,const InterfaceGeom& iface,std::mt19937_64& rng) {
    CDT cdt;
    const double h=cfg.target_edge_length;
    // Case 6 can approach a thin strip above edge 0 when the right-edge
    // crossing approaches corner 1. Resolve that strip using spacing tied
    // to its physical thickness instead of relaxing the quality threshold.
    double thin_h=h;
    if(cid==6){
        const double thickness=2.0*std::clamp(iface.pb,0.0,1.0);
        if(thickness>0.0)thin_h=std::min(h,30.0*thickness);
    }
    std::array<std::set<double>,4> bt;
    for(int e=0;e<4;++e){
        const double eh=(cid==6 && e==0)?thin_h:h;
        const int n=std::max(1,static_cast<int>(std::ceil(2.0/eh)));
        for(int i=0;i<=n;++i)bt[static_cast<std::size_t>(e)].insert(static_cast<double>(i)/n);
    }
    auto grade_boundary=[&](P2 p,int e){
        if(e<0)return;
        const double t=edge_parameter(e,p);
        bt[static_cast<std::size_t>(e)].insert(std::clamp(t,0.0,1.0));
        const double d=std::min(t,1.0-t);
        if(!(d>0.0))return;
        const bool near_start=t<0.5;
        const int adjacent=near_start ? (e+3)%4 : (e+1)%4;
        double s=d;
        const double target=std::min(0.25,h/2.0);
        while(s<target){
            if(t-s>0.0)bt[static_cast<std::size_t>(e)].insert(t-s);
            if(t+s<1.0)bt[static_cast<std::size_t>(e)].insert(t+s);
            // Grade the other square edge incident to the same near corner.
            const double ta=near_start ? 1.0-s : s;
            if(ta>0.0&&ta<1.0)bt[static_cast<std::size_t>(adjacent)].insert(ta);
            s*=2.0;
        }
    };
    grade_boundary(iface.a,iface.edge_a);
    grade_boundary(iface.b,iface.edge_b);

    std::vector<std::vector<VH>> boundary(4);
    for(int e=0;e<4;++e){
        for(double t:bt[static_cast<std::size_t>(e)]){
            const auto p=edge_point(e,t);
            boundary[static_cast<std::size_t>(e)].push_back(cdt.insert(Point(p.x,p.y)));
        }
        add_constraint_chain(cdt,boundary[static_cast<std::size_t>(e)]);
    }

    std::set<double> it;
    const double L=std::sqrt(dist2(iface.a,iface.b));
    const double ih=(cid==6)?thin_h:h;
    const int ni=std::max(1,static_cast<int>(std::ceil(L/ih)));
    for(int i=0;i<=ni;++i)it.insert(static_cast<double>(i)/ni);
    auto grade_interface=[&](P2 p,int e,bool from_a){
        if(e<0||!(L>0.0))return;
        const double t=edge_parameter(e,p);
        const double physical=2.0*std::min(t,1.0-t);
        if(!(physical>0.0))return;
        double tau=physical/L;
        while(tau<1.0){
            if(from_a)it.insert(std::min(1.0,tau)); else it.insert(std::max(0.0,1.0-tau));
            if(tau*L>=h)break;
            tau*=2.0;
        }
    };
    grade_interface(iface.a,iface.edge_a,true);
    grade_interface(iface.b,iface.edge_b,false);
    std::vector<VH> ichain;
    for(double t:it){
        P2 p{iface.a.x+t*(iface.b.x-iface.a.x),iface.a.y+t*(iface.b.y-iface.a.y)};
        ichain.push_back(cdt.insert(Point(p.x,p.y)));
    }
    add_constraint_chain(cdt,ichain);

    int nx=std::max(2,static_cast<int>(std::ceil(2.0/h)));
    if(cid==6 && thin_h<h)nx=std::min(nx,5);
    const double step=2.0/nx;
    std::uniform_real_distribution<double> jitter(-0.22,0.22);
    for (int j=1;j<nx;++j) for (int i=1;i<nx;++i) {
        double x=-1.0+i*step+step*jitter(rng);
        double y=-1.0+j*step+step*jitter(rng);
        x=std::clamp(x,-1.0+0.05*step,1.0-0.05*step);
        y=std::clamp(y,-1.0+0.05*step,1.0-0.05*step);
        P2 p{x,y};
        const double d=std::abs(phi_raw(iface,p));
        if (d<0.10*step) continue;
        cdt.insert(Point(x,y));
    }

    // Bounded quality refinement. Poor finite triangles receive a centroid
    // insertion when that point is safely away from the interface/boundary.
    const std::size_t primary_budget=std::max<std::size_t>(4,cfg.max_nodes/3);
    for(int refine=0;refine<64 && cdt.number_of_vertices()<primary_budget;++refine){
        double worst=1.0;
        std::optional<P2> candidate;
        for(auto f=cdt.finite_faces_begin();f!=cdt.finite_faces_end();++f){
            P2 a{CGAL::to_double(f->vertex(0)->point().x()),CGAL::to_double(f->vertex(0)->point().y())};
            P2 b{CGAL::to_double(f->vertex(1)->point().x()),CGAL::to_double(f->vertex(1)->point().y())};
            P2 c{CGAL::to_double(f->vertex(2)->point().x()),CGAL::to_double(f->vertex(2)->point().y())};
            if(cross(a,b,c)<=0)std::swap(b,c);
            const double q=triangle_quality(a,b,c);
            if(q>=cfg.min_triangle_quality*1.05 || q>=worst)continue;
            P2 p{(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
            if(p.x<=-1.0+1e-8||p.x>=1.0-1e-8||p.y<=-1.0+1e-8||p.y>=1.0-1e-8)continue;
            if(std::abs(phi_raw(iface,p))<0.08*step)continue;
            worst=q;candidate=p;
        }
        if(!candidate)break;
        cdt.insert(Point(candidate->x,candidate->y));
    }

    Mesh mesh; mesh.iface=iface;
    std::map<const void*,int> ids;
    auto idof=[&](VH vh)->int{
        const void* key=static_cast<const void*>(&*vh);
        auto it=ids.find(key); if(it!=ids.end()) return it->second;
        P2 p{CGAL::to_double(vh->point().x()),CGAL::to_double(vh->point().y())};
        if (!inside_square(p)) throw std::runtime_error("triangulation vertex outside square");
        Node n=make_node(p,iface,cfg,true); n.id=static_cast<int>(mesh.nodes.size());
        mesh.nodes.push_back(n); ids[key]=n.id; return n.id;
    };

    for (auto f=cdt.finite_faces_begin();f!=cdt.finite_faces_end();++f) {
        P2 p[3]; for(int k=0;k<3;++k) p[k]={CGAL::to_double(f->vertex(k)->point().x()),CGAL::to_double(f->vertex(k)->point().y())};
        P2 c{(p[0].x+p[1].x+p[2].x)/3.0,(p[0].y+p[1].y+p[2].y)/3.0};
        if (!inside_square(c)) continue;
        Tri t; for(int k=0;k<3;++k)t.v[static_cast<std::size_t>(k)]=idof(f->vertex(k));
        if (cross(mesh.nodes[t.v[0]].p,mesh.nodes[t.v[1]].p,mesh.nodes[t.v[2]].p)<0) std::swap(t.v[1],t.v[2]);
        const double ph=phi_raw(iface,c);
        t.phase=ph<0?-1:1;
        t.quality=triangle_quality(mesh.nodes[t.v[0]].p,mesh.nodes[t.v[1]].p,mesh.nodes[t.v[2]].p);
        mesh.tris.push_back(t);
    }
    if (mesh.nodes.size()>cfg.max_nodes) throw std::runtime_error("primary node count exceeds max-nodes");
    if (mesh.tris.empty()) throw std::runtime_error("empty triangulation");

    // Ordered constrained edge lists are reconstructed geometrically from primary nodes.
    for (int e=0;e<4;++e) {
        std::vector<std::pair<double,int>> v;
        for (const auto& n:mesh.nodes) if (boundary_edge(n.p)==e) v.push_back({edge_parameter(e,n.p),n.id});
        std::sort(v.begin(),v.end());
        for(std::size_t i=1;i<v.size();++i) mesh.boundary_edges.push_back({v[i-1].second,v[i].second});
    }
    std::vector<std::pair<double,int>> iv;
    for (const auto& n:mesh.nodes) if (std::abs(phi_raw(iface,n.p))<=5e-10) {
        const double t=((n.p.x-iface.a.x)*(iface.b.x-iface.a.x)+(n.p.y-iface.a.y)*(iface.b.y-iface.a.y))/dist2(iface.a,iface.b);
        if(t>=-kGeomTol&&t<=1+kGeomTol) iv.push_back({t,n.id});
    }
    std::sort(iv.begin(),iv.end());
    for(std::size_t i=1;i<iv.size();++i) mesh.interface_edges.push_back({iv[i-1].second,iv[i].second});
    return mesh;
}

std::array<int,4> quad_from_pair(const Mesh& m,const Tri& a,const Tri& b) {
    std::set<int> u(a.v.begin(),a.v.end()); u.insert(b.v.begin(),b.v.end());
    if(u.size()!=4) return {-1,-1,-1,-1};
    P2 c{};
    for(int id:u){c.x+=m.nodes[id].p.x;c.y+=m.nodes[id].p.y;} c.x/=4;c.y/=4;
    std::vector<int> ids(u.begin(),u.end());
    std::sort(ids.begin(),ids.end(),[&](int i,int j){
        return std::atan2(m.nodes[i].p.y-c.y,m.nodes[i].p.x-c.x) < std::atan2(m.nodes[j].p.y-c.y,m.nodes[j].p.x-c.x);
    });
    std::array<int,4> q{ids[0],ids[1],ids[2],ids[3]};
    std::array<P2,4> p{m.nodes[q[0]].p,m.nodes[q[1]].p,m.nodes[q[2]].p,m.nodes[q[3]].p};
    if(area_poly(std::vector<P2>(p.begin(),p.end()))<0){std::reverse(q.begin(),q.end());}
    return q;
}
bool convex_quad(const Mesh& m,const std::array<int,4>& q) {
    if(q[0]<0)return false;
    double s=0;
    for(int i=0;i<4;++i){
        const double c=cross(m.nodes[q[i]].p,m.nodes[q[(i+1)%4]].p,m.nodes[q[(i+2)%4]].p);
        if(c<=kGeomTol)return false; s+=c;
    }
    return s>0;
}
bool constrained_shared_edge(const Mesh& m,int a,int b) {
    const auto& x=m.nodes[a]; const auto& y=m.nodes[b];
    if (std::abs(phi_raw(m.iface,x.p))<=5e-10 && std::abs(phi_raw(m.iface,y.p))<=5e-10) return true;
    const int ex=boundary_edge(x.p), ey=boundary_edge(y.p);
    return ex>=0&&ex==ey;
}

void pair_triangles(Mesh& m,const Config& cfg) {
    struct EInfo{int a{-1},b{-1};};
    std::map<std::array<int,2>,EInfo> emap;
    for(int ti=0;ti<static_cast<int>(m.tris.size());++ti){
        for(int e=0;e<3;++e){
            auto k=edge_key(m.tris[ti].v[e],m.tris[ti].v[(e+1)%3]);
            auto& x=emap[k]; if(x.a<0)x.a=ti;else x.b=ti;
        }
    }
    using G=boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,boost::property<boost::edge_weight_t,long long>>;
    G g(m.tris.size());
    constexpr long long Q=100000;
    const long long max_pairs=static_cast<long long>(m.tris.size()/2);
    if(max_pairs>(std::numeric_limits<long long>::max()/(Q+1))-2) throw std::overflow_error("matching weight overflow");
    const long long C=(max_pairs+1)*(Q+1);
    struct Cand{int a,b;std::array<int,4> q;double quality;};
    std::map<std::array<int,2>,Cand> cand;
    for(const auto& [edge,ei]:emap){
        if(ei.a<0||ei.b<0)continue;
        const auto& ta=m.tris[ei.a]; const auto& tb=m.tris[ei.b];
        if(ta.phase!=tb.phase||constrained_shared_edge(m,edge[0],edge[1]))continue;
        auto q=quad_from_pair(m,ta,tb);
        if(!convex_quad(m,q))continue;
        std::array<P2,4> qp{m.nodes[q[0]].p,m.nodes[q[1]].p,m.nodes[q[2]].p,m.nodes[q[3]].p};
        const double qual=quad_quality(qp);
        if(qual+1e-14<cfg.min_quad_quality)continue;
        const long long bonus=static_cast<long long>(std::llround(std::clamp(qual,0.0,1.0)*Q));
        if(C>std::numeric_limits<long long>::max()-bonus)throw std::overflow_error("matching weight overflow");
        boost::add_edge(static_cast<std::size_t>(ei.a),static_cast<std::size_t>(ei.b),C+bonus,g);
        cand[edge_key(ei.a,ei.b)]={ei.a,ei.b,q,qual};
    }
    std::vector<boost::graph_traits<G>::vertex_descriptor> mate(boost::num_vertices(g));
    boost::maximum_weighted_matching(g,boost::make_iterator_property_map(mate.begin(),get(boost::vertex_index,g)));
    std::vector<bool> used(m.tris.size(),false);
    std::vector<Quad> quads; std::vector<Tri> keep;
    const auto nullv=boost::graph_traits<G>::null_vertex();
    for(std::size_t i=0;i<mate.size();++i){
        const auto j=mate[i]; if(j==nullv||i>=j)continue;
        auto it=cand.find(edge_key(static_cast<int>(i),static_cast<int>(j)));
        if(it==cand.end())continue;
        used[i]=used[j]=true;
        Quad q; q.v=it->second.q; q.phase=m.tris[i].phase; q.quality=it->second.quality; quads.push_back(q);
    }
    for(std::size_t i=0;i<m.tris.size();++i)if(!used[i])keep.push_back(m.tris[i]);
    m.tris=std::move(keep); m.quads=std::move(quads);
}
double local_min_quality(const Mesh& m,int nid) {
    double q=std::numeric_limits<double>::infinity(); bool any=false;
    for(const auto& t:m.tris)if(std::find(t.v.begin(),t.v.end(),nid)!=t.v.end()){
        q=std::min(q,triangle_quality(m.nodes[t.v[0]].p,m.nodes[t.v[1]].p,m.nodes[t.v[2]].p));any=true;
    }
    for(const auto& e:m.quads)if(std::find(e.v.begin(),e.v.end(),nid)!=e.v.end()){
        std::array<P2,4> p{m.nodes[e.v[0]].p,m.nodes[e.v[1]].p,m.nodes[e.v[2]].p,m.nodes[e.v[3]].p};
        q=std::min(q,quad_quality(p));any=true;
    }
    return any?q:1.0;
}
std::vector<int> neighbors(const Mesh& m,int nid){
    std::set<int>s;
    for(const auto&t:m.tris)if(std::find(t.v.begin(),t.v.end(),nid)!=t.v.end())for(int x:t.v)if(x!=nid)s.insert(x);
    for(const auto&q:m.quads)if(std::find(q.v.begin(),q.v.end(),nid)!=q.v.end())for(int x:q.v)if(x!=nid)s.insert(x);
    return {s.begin(),s.end()};
}
P2 project_constraint(const Node& n,P2 p,const InterfaceGeom& g){
    if(n.constraint=="interior") return p;
    if(n.constraint=="boundary"){
        const auto a=corners[static_cast<std::size_t>(n.parent_edge)],b=corners[static_cast<std::size_t>((n.parent_edge+1)%4)];
        const double dx=b.x-a.x,dy=b.y-a.y;
        double t=((p.x-a.x)*dx+(p.y-a.y)*dy)/(dx*dx+dy*dy); t=std::clamp(t,0.0,1.0);
        return {a.x+t*dx,a.y+t*dy};
    }
    if(n.constraint=="interface"){
        const double dx=g.b.x-g.a.x,dy=g.b.y-g.a.y;
        double t=((p.x-g.a.x)*dx+(p.y-g.a.y)*dy)/(dx*dx+dy*dy); t=std::clamp(t,0.0,1.0);
        return {g.a.x+t*dx,g.a.y+t*dy};
    }
    return n.p;
}
bool ordering_ok(const Mesh&m,int nid,P2 p){
    const auto& n=m.nodes[nid];
    if(n.constraint!="boundary"&&n.constraint!="interface")return true;
    auto param=[&](const Node& x,P2 pt){
        if(x.constraint=="boundary")return edge_parameter(x.parent_edge,pt);
        const double dx=m.iface.b.x-m.iface.a.x,dy=m.iface.b.y-m.iface.a.y;
        return ((pt.x-m.iface.a.x)*dx+(pt.y-m.iface.a.y)*dy)/(dx*dx+dy*dy);
    };
    const double t=param(n,p); double lo=0.0,hi=1.0;
    for(const auto& o:m.nodes){
        if(o.id==nid||!o.primary)continue;
        bool same=false;
        if(n.constraint=="boundary")same=o.constraint=="boundary"&&o.parent_edge==n.parent_edge;
        else same=o.constraint=="interface"&&o.parent_segment==n.parent_segment;
        if(!same)continue;
        const double ot=param(o,o.p);
        if(ot<n.parameter)lo=std::max(lo,ot); else if(ot>n.parameter)hi=std::min(hi,ot);
    }
    return t>lo+1e-10&&t<hi-1e-10;
}
void refresh_quality(Mesh&m){
    for(auto&t:m.tris)t.quality=triangle_quality(m.nodes[t.v[0]].p,m.nodes[t.v[1]].p,m.nodes[t.v[2]].p);
    for(auto&q:m.quads){
        std::array<P2,4> p{m.nodes[q.v[0]].p,m.nodes[q.v[1]].p,m.nodes[q.v[2]].p,m.nodes[q.v[3]].p};
        q.quality=quad_quality(p);
    }
}
void smooth(Mesh&m,const Config&cfg){
    for(int pass=0;pass<cfg.max_smoothing_passes;++pass){
        bool changed=false;
        for(int nid=0;nid<static_cast<int>(m.nodes.size());++nid){
            auto& n=m.nodes[nid]; if(!n.primary||n.fixed||n.constraint=="corner")continue;
            const auto nb=neighbors(m,nid); if(nb.size()<2)continue;
            P2 cen{};for(int j:nb){cen.x+=m.nodes[j].p.x;cen.y+=m.nodes[j].p.y;}cen.x/=nb.size();cen.y/=nb.size();
            cen=project_constraint(n,cen,m.iface);
            const P2 old=n.p; const double before=local_min_quality(m,nid);
            double alpha=1.0;
            for(int ls=0;ls<12;++ls){
                P2 c{old.x+alpha*(cen.x-old.x),old.y+alpha*(cen.y-old.y)};
                c=project_constraint(n,c,m.iface);
                if(!inside_square(c)||!ordering_ok(m,nid,c)){alpha*=0.5;continue;}
                n.p=c;
                const double after=local_min_quality(m,nid);
                if(std::isfinite(after)&&after>before+1e-12){changed=true;break;}
                n.p=old; alpha*=0.5;
            }
            if(dist2(n.p,old)<1e-30)n.p=old;
        }
        if(!changed)break;
    }
    for(auto&n:m.nodes){
        if(n.constraint=="interface")n.phi=0.0;
        else n.phi=phi_raw(m.iface,n.p);
        n.sign=classify_phi(n.phi,cfg.phi_zero_tol);
        if(n.constraint=="boundary")n.parameter=edge_parameter(n.parent_edge,n.p);
        else if(n.constraint=="interface"){
            const double dx=m.iface.b.x-m.iface.a.x,dy=m.iface.b.y-m.iface.a.y;
            n.parameter=((n.p.x-m.iface.a.x)*dx+(n.p.y-m.iface.a.y)*dy)/(dx*dx+dy*dy);
        }
    }
    refresh_quality(m);
}
std::vector<P2> clip_halfplane(std::vector<P2> poly,const InterfaceGeom&g,bool negative){
    std::vector<P2> out;
    auto val=[&](P2 p){double v=phi_raw(g,p);return negative?v:-v;};
    for(std::size_t i=0;i<poly.size();++i){
        P2 a=poly[i],b=poly[(i+1)%poly.size()]; double fa=val(a),fb=val(b);
        bool ia=fa<=0,ib=fb<=0;
        if(ia)out.push_back(a);
        if(ia!=ib){
            const double t=fa/(fa-fb);
            out.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});
        }
    }
    return out;
}
void update_metrics(Mesh&m){
    m.tri_min=m.tris.empty()?1.0:std::numeric_limits<double>::infinity();
    m.quad_min=m.quads.empty()?1.0:std::numeric_limits<double>::infinity();
    double ts=0,qs=0,qa=0,ta=0;
    m.phase_area={0,0};
    for(const auto&t:m.tris){
        m.tri_min=std::min(m.tri_min,t.quality);ts+=t.quality;
        const double a=0.5*cross(m.nodes[t.v[0]].p,m.nodes[t.v[1]].p,m.nodes[t.v[2]].p);ta+=a;
        m.phase_area[t.phase<0?0:1]+=a;
    }
    for(const auto&q:m.quads){
        m.quad_min=std::min(m.quad_min,q.quality);qs+=q.quality;
        std::array<P2,4> p{m.nodes[q.v[0]].p,m.nodes[q.v[1]].p,m.nodes[q.v[2]].p,m.nodes[q.v[3]].p};
        const double a=quad_area(p);qa+=a;m.phase_area[q.phase<0?0:1]+=a;
    }
    m.tri_mean=m.tris.empty()?1.0:ts/m.tris.size();
    m.quad_mean=m.quads.empty()?1.0:qs/m.quads.size();
    const double ne=static_cast<double>(m.tris.size()+m.quads.size());
    m.quad_count_fraction=ne?m.quads.size()/ne:0;
    m.quad_area_fraction=(ta+qa)>0?qa/(ta+qa):0;
}
void enrich(Mesh&m,const Config&cfg){
    std::map<std::array<int,2>,int> mids;
    auto midpoint=[&](int a,int b){
        auto k=edge_key(a,b);auto it=mids.find(k);if(it!=mids.end())return it->second;
        P2 p{0.5*(m.nodes[a].p.x+m.nodes[b].p.x),0.5*(m.nodes[a].p.y+m.nodes[b].p.y)};
        Node n=make_node(p,m.iface,cfg,false);n.id=static_cast<int>(m.nodes.size());
        if(std::abs(m.nodes[a].phi)<=cfg.phi_zero_tol&&std::abs(m.nodes[b].phi)<=cfg.phi_zero_tol){n.constraint="interface";n.phi=0;n.sign=0;n.parent_segment=0;}
        m.nodes.push_back(n);mids[k]=n.id;return n.id;
    };
    for(auto&t:m.tris){
        t.p2={t.v[0],t.v[1],t.v[2],midpoint(t.v[0],t.v[1]),midpoint(t.v[1],t.v[2]),midpoint(t.v[2],t.v[0])};
    }
    for(auto&q:m.quads){
        q.q2[0]=q.v[0];q.q2[1]=q.v[1];q.q2[2]=q.v[2];q.q2[3]=q.v[3];
        q.q2[4]=midpoint(q.v[0],q.v[1]);q.q2[5]=midpoint(q.v[1],q.v[2]);q.q2[6]=midpoint(q.v[2],q.v[3]);q.q2[7]=midpoint(q.v[3],q.v[0]);
        P2 p{};for(int id:q.v){p.x+=m.nodes[id].p.x;p.y+=m.nodes[id].p.y;}p.x/=4;p.y/=4;
        Node n=make_node(p,m.iface,cfg,false);n.id=static_cast<int>(m.nodes.size());m.nodes.push_back(n);q.q2[8]=n.id;
    }
    if(m.nodes.size()>cfg.max_nodes)throw std::runtime_error("enriched node count exceeds max-nodes");

    std::map<std::pair<int,int>,int> pid;
    auto pget=[&](int geom,int phase){
        auto k=std::make_pair(geom,phase);auto it=pid.find(k);if(it!=pid.end())return it->second;
        int id=static_cast<int>(m.pressures.size());m.pressures.push_back({id,geom,phase});pid[k]=id;return id;
    };
    for(auto&t:m.tris)for(int i=0;i<3;++i)t.pressure[i]=pget(t.v[i],t.phase);
    for(auto&q:m.quads)for(int i=0;i<4;++i)q.pressure[i]=pget(q.v[i],q.phase);
}
bool cell_crosses_interface(const Mesh&m,const std::vector<int>&v,int phase,const Config&cfg){
    for(int id:v){
        const double ph=phi_raw(m.iface,m.nodes[id].p);
        const int s=classify_phi(ph,cfg.phi_zero_tol*2);
        if(s!=0&&s!=phase)return true;
    }
    return false;
}
void validate_mesh(Mesh&m,const Config&cfg,int cid){
    if(m.nodes.empty())throw std::runtime_error("no nodes");
    const auto want=expected_signs(cid);
    for(int i=0;i<4;++i){
        const double ph=phi_raw(m.iface,corners[static_cast<std::size_t>(i)]);
        m.corner_phi[i]=std::abs(ph)<=cfg.phi_zero_tol?0.0:ph;
        m.corner_sign[i]=classify_phi(m.corner_phi[i],cfg.phi_zero_tol);
        if(m.corner_sign[i]!=want[i])throw std::runtime_error("canonical corner sign mismatch");
    }
    for(const auto&n:m.nodes){
        if(!std::isfinite(n.p.x)||!std::isfinite(n.p.y)||!std::isfinite(n.phi))throw std::runtime_error("non-finite node");
        if(!inside_square(n.p))throw std::runtime_error("node outside square");
        if(n.constraint=="boundary"&&boundary_edge(n.p)!=n.parent_edge)throw std::runtime_error("boundary constraint violation");
        if(n.constraint=="interface"&&std::abs(phi_raw(m.iface,n.p))>1e-8)throw std::runtime_error("interface constraint violation");
    }
    std::set<std::vector<int>> cells;
    std::map<std::array<int,2>,int> incidence;
    double area=0;
    for(const auto&t:m.tris){
        const double q=triangle_quality(m.nodes[t.v[0]].p,m.nodes[t.v[1]].p,m.nodes[t.v[2]].p);
        if(q+1e-12<cfg.min_triangle_quality)throw std::runtime_error("triangle quality threshold violated");
        std::vector<int> vv(t.v.begin(),t.v.end()); if(cell_crosses_interface(m,vv,t.phase,cfg))throw std::runtime_error("triangle crosses interface");
        auto sig=vv;std::sort(sig.begin(),sig.end());if(!cells.insert(sig).second)throw std::runtime_error("duplicate cell");
        const double a=0.5*cross(m.nodes[t.v[0]].p,m.nodes[t.v[1]].p,m.nodes[t.v[2]].p);if(!(a>0))throw std::runtime_error("inverted triangle");area+=a;
        for(int e=0;e<3;++e)++incidence[edge_key(t.v[e],t.v[(e+1)%3])];
    }
    for(const auto&q:m.quads){
        std::array<P2,4> p{m.nodes[q.v[0]].p,m.nodes[q.v[1]].p,m.nodes[q.v[2]].p,m.nodes[q.v[3]].p};
        const double qual=quad_quality(p);if(qual+1e-12<cfg.min_quad_quality)throw std::runtime_error("quad quality threshold violated");
        if(!convex_quad(m,q.v))throw std::runtime_error("nonconvex quad");
        std::vector<int> vv(q.v.begin(),q.v.end());if(cell_crosses_interface(m,vv,q.phase,cfg))throw std::runtime_error("quad crosses interface");
        auto sig=vv;std::sort(sig.begin(),sig.end());if(!cells.insert(sig).second)throw std::runtime_error("duplicate cell");
        area+=quad_area(p);for(int e=0;e<4;++e)++incidence[edge_key(q.v[e],q.v[(e+1)%4])];
    }
    if(std::abs(area-4.0)>kAreaTol)throw std::runtime_error("mesh has holes/overlaps: total area mismatch");
    for(const auto&[e,c]:incidence)if(c<1||c>2)throw std::runtime_error("invalid shared-edge incidence");
    const auto neg=clip_halfplane({corners.begin(),corners.end()},m.iface,true);
    const double an=std::abs(area_poly(neg)),ap=4.0-an;
    if(std::abs(m.phase_area[0]-an)>5e-7||std::abs(m.phase_area[1]-ap)>5e-7)throw std::runtime_error("phase areas disagree with independent clipping");

    // Interface geometric vertices used on both sides must have distinct pressure records.
    std::map<int,std::set<int>> phase_by_geom;
    for(const auto&p:m.pressures)phase_by_geom[p.geom].insert(p.phase);
    for(const auto&n:m.nodes)if(n.primary&&n.constraint=="interface"){
        std::set<int> incident;
        for(const auto&t:m.tris)if(std::find(t.v.begin(),t.v.end(),n.id)!=t.v.end())incident.insert(t.phase);
        for(const auto&q:m.quads)if(std::find(q.v.begin(),q.v.end(),n.id)!=q.v.end())incident.insert(q.phase);
        if(incident.size()==2&&phase_by_geom[n.id].size()!=2)throw std::runtime_error("pressure-side separation missing");
    }
}
std::string serialize(const Mesh&m,const Config&cfg,int cid,int tid,std::uint64_t seed){
    std::ostringstream o;o<<std::setprecision(17);
    o<<"PVMLS_TEMPLATE_DATA "<<kSchemaVersion<<"\n";
    o<<"META\nCASE "<<cid<<"\nTEMPLATE "<<tid<<"\nSEED "<<seed<<"\n";
    o<<"CONFIG phi_zero_tol "<<cfg.phi_zero_tol<<" min_edge_fraction "<<cfg.min_edge_fraction
     <<" target_edge_length "<<cfg.target_edge_length<<" min_triangle_quality "<<cfg.min_triangle_quality
     <<" min_quad_quality "<<cfg.min_quad_quality<<" max_nodes "<<cfg.max_nodes
     <<" max_attempts_per_template "<<cfg.max_attempts_per_template<<" max_smoothing_passes "<<cfg.max_smoothing_passes
     <<" image_size "<<cfg.image_size<<"\n";
    o<<"CORNERS 4\n";
    for(int i=0;i<4;++i)o<<i<<" "<<corners[i].x<<" "<<corners[i].y<<" "<<m.corner_phi[i]<<" "<<m.corner_sign[i]<<"\n";
    o<<"INTERFACE\n";
    o<<"INTERFACE_ENDPOINT 0 "<<m.iface.a.x<<" "<<m.iface.a.y<<" edge "<<m.iface.edge_a<<" parameter "<<m.iface.pa<<"\n";
    o<<"INTERFACE_ENDPOINT 1 "<<m.iface.b.x<<" "<<m.iface.b.y<<" edge "<<m.iface.edge_b<<" parameter "<<m.iface.pb<<"\n";
    o<<"INTERFACE_SIGN_FACTOR "<<m.iface.sign_factor<<"\n";
    o<<"NODES "<<m.nodes.size()<<"\n";
    o<<"# id ksi eta phi sign constraint parent_edge parent_segment parameter primary\n";
    for(const auto&n:m.nodes)o<<n.id<<" "<<n.p.x<<" "<<n.p.y<<" "<<n.phi<<" "<<n.sign<<" "<<n.constraint<<" "<<n.parent_edge<<" "<<n.parent_segment<<" "<<n.parameter<<" "<<(n.primary?1:0)<<"\n";
    o<<"PRESSURE_RECORDS "<<m.pressures.size()<<"\n";
    for(const auto&p:m.pressures)o<<p.id<<" "<<p.geom<<" "<<p.phase<<"\n";
    o<<"TRIANGLES "<<m.tris.size()<<"\n";
    for(std::size_t i=0;i<m.tris.size();++i){const auto&t=m.tris[i];o<<i<<" phase "<<t.phase<<" quality "<<t.quality<<" P1";for(int x:t.v)o<<" "<<x;o<<" P2";for(int x:t.p2)o<<" "<<x;o<<"\n";}
    o<<"QUADS "<<m.quads.size()<<"\n";
    for(std::size_t i=0;i<m.quads.size();++i){const auto&q=m.quads[i];o<<i<<" phase "<<q.phase<<" quality "<<q.quality<<" Q1";for(int x:q.v)o<<" "<<x;o<<" Q2";for(int x:q.q2)o<<" "<<x;o<<"\n";}
    o<<"PRESSURE_CONNECTIVITY\n";
    for(std::size_t i=0;i<m.tris.size();++i){o<<"T "<<i;for(int x:m.tris[i].pressure)o<<" "<<x;o<<"\n";}
    for(std::size_t i=0;i<m.quads.size();++i){o<<"Q "<<i;for(int x:m.quads[i].pressure)o<<" "<<x;o<<"\n";}
    o<<"BOUNDARY_EDGES "<<m.boundary_edges.size()<<"\n";for(const auto&e:m.boundary_edges)o<<e[0]<<" "<<e[1]<<"\n";
    o<<"INTERFACE_EDGES "<<m.interface_edges.size()<<"\n";for(const auto&e:m.interface_edges)o<<e[0]<<" "<<e[1]<<"\n";
    o<<"METRICS tri_min "<<m.tri_min<<" tri_mean "<<m.tri_mean<<" quad_min "<<m.quad_min<<" quad_mean "<<m.quad_mean
     <<" quad_count_fraction "<<m.quad_count_fraction<<" quad_area_fraction "<<m.quad_area_fraction
     <<" phase_negative_area "<<m.phase_area[0]<<" phase_positive_area "<<m.phase_area[1]<<"\n";
    o<<"VALIDATION OK\nEND\n";return o.str();
}
void validate_serialized_roundtrip(const std::string&dat,const Mesh&m){
    std::istringstream in(dat);
    std::string line;
    std::size_t count=0;
    bool nodes=false;
    while(std::getline(in,line)){
        if(line.rfind("NODES ",0)==0){nodes=true;continue;}
        if(!nodes)continue;
        if(line.empty()||line[0]=='#')continue;
        if(line.rfind("PRESSURE_RECORDS ",0)==0)break;
        std::istringstream row(line);
        int id=-1, sign=0, pe=-1, ps=-1, primary=0;
        double x=0,y=0,phi=0,param=0;
        std::string constraint;
        if(!(row>>id>>x>>y>>phi>>sign>>constraint>>pe>>ps>>param>>primary))
            throw std::runtime_error("serialized node parse failed");
        if(id<0||static_cast<std::size_t>(id)>=m.nodes.size())throw std::runtime_error("serialized node id out of range");
        const auto& n=m.nodes[static_cast<std::size_t>(id)];
        if(std::bit_cast<std::uint64_t>(x)!=std::bit_cast<std::uint64_t>(n.p.x) ||
           std::bit_cast<std::uint64_t>(y)!=std::bit_cast<std::uint64_t>(n.p.y) ||
           std::bit_cast<std::uint64_t>(phi)!=std::bit_cast<std::uint64_t>(n.phi))
            throw std::runtime_error("serialized floating-point round trip failed");
        ++count;
    }
    if(count!=m.nodes.size())throw std::runtime_error("serialized node count mismatch");
}
std::string signature_of(const Mesh&m){
    std::ostringstream o;o<<std::setprecision(14);
    o<<m.iface.a.x<<","<<m.iface.a.y<<","<<m.iface.b.x<<","<<m.iface.b.y<<";";
    for(const auto&n:m.nodes)if(n.primary)o<<n.p.x<<","<<n.p.y<<";";
    for(const auto&t:m.tris){for(int x:t.v)o<<x<<",";o<<"T";}
    for(const auto&q:m.quads){for(int x:q.v)o<<x<<",";o<<"Q";}
    return o.str();
}
Candidate build_candidate(const Config&cfg,int cid,int tid,int attempt){
    std::mt19937_64 rng(mix_seed(cfg.seed,cid,tid,attempt));
    const auto iface=sample_interface(cfg,cid,tid,attempt,rng);
    Mesh m=triangulate_primary(cfg,cid,iface,rng);
    pair_triangles(m,cfg);
    smooth(m,cfg);
    refresh_quality(m);update_metrics(m);
    enrich(m,cfg);update_metrics(m);
    validate_mesh(m,cfg,cid);
    return {std::move(m),{}};
}

// Tiny 5x7 glyph set used only for preview identifiers/legend.
std::array<unsigned char,7> glyph(char c){
    switch(c){
    case 'A':return{14,17,17,31,17,17,17}; case 'B':return{30,17,17,30,17,17,30}; case 'C':return{14,17,16,16,16,17,14};
    case 'E':return{31,16,16,30,16,16,31}; case 'G':return{14,17,16,23,17,17,15};
    case 'I':return{31,4,4,4,4,4,31}; case 'L':return{16,16,16,16,16,16,31};
    case 'N':return{17,25,21,19,17,17,17}; case 'P':return{30,17,17,30,16,16,16};
    case 'O':return{14,17,17,17,17,17,14}; case 'Q':return{14,17,17,17,21,18,13};
    case 'R':return{30,17,17,30,20,18,17}; case 'S':return{15,16,16,14,1,1,30};
    case 'T':return{31,4,4,4,4,4,4}; case 'U':return{17,17,17,17,17,17,14};
    case 'D':return{30,17,17,17,17,17,30}; case 'M':return{17,27,21,17,17,17,17};
    case '0':return{14,17,19,21,25,17,14}; case '1':return{4,12,4,4,4,4,14};
    case '2':return{14,17,1,2,4,8,31}; case '3':return{30,1,1,14,1,1,30};
    case '4':return{2,6,10,18,31,2,2}; case '5':return{31,16,16,30,1,1,30};
    case '6':return{14,16,16,30,17,17,14}; case '7':return{31,1,2,4,8,8,8};
    case '8':return{14,17,17,14,17,17,14}; case '9':return{14,17,17,15,1,1,14};
    case '-':return{0,0,0,31,0,0,0}; default:return{0,0,0,0,0,0,0};
    }
}
struct Image{
    int w,h;std::vector<unsigned char>px;
    Image(int W,int H):w(W),h(H),px(static_cast<std::size_t>(W*H*3),245){}
    void set(int x,int y,unsigned char r,unsigned char g,unsigned char b){if(x<0||y<0||x>=w||y>=h)return;auto i=static_cast<std::size_t>((y*w+x)*3);px[i]=r;px[i+1]=g;px[i+2]=b;}
    void line(int x0,int y0,int x1,int y1,unsigned char r,unsigned char g,unsigned char b,int thick=1){
        int dx=std::abs(x1-x0),sx=x0<x1?1:-1,dy=-std::abs(y1-y0),sy=y0<y1?1:-1,err=dx+dy;
        for(;;){for(int oy=-thick/2;oy<=thick/2;++oy)for(int ox=-thick/2;ox<=thick/2;++ox)set(x0+ox,y0+oy,r,g,b);if(x0==x1&&y0==y1)break;int e2=2*err;if(e2>=dy){err+=dy;x0+=sx;}if(e2<=dx){err+=dx;y0+=sy;}}
    }
    void text(int x,int y,const std::string&s,int scale=2){for(char c:s){auto g=glyph(c);for(int yy=0;yy<7;++yy)for(int xx=0;xx<5;++xx)if(g[yy]&(1u<<(4-xx)))for(int a=0;a<scale;++a)for(int b=0;b<scale;++b)set(x+xx*scale+a,y+yy*scale+b,20,20,20);x+=6*scale;}}
};
void fill_poly(Image&im,const std::vector<std::pair<int,int>>&p,unsigned char r,unsigned char g,unsigned char b){
    if(p.size()<3)return;int ymin=im.h-1,ymax=0;for(auto[x,y]:p){ymin=std::min(ymin,y);ymax=std::max(ymax,y);}
    for(int y=std::max(0,ymin);y<=std::min(im.h-1,ymax);++y){
        std::vector<int> xs;
        for(std::size_t i=0;i<p.size();++i){auto[x1,y1]=p[i];auto[x2,y2]=p[(i+1)%p.size()];if((y1<=y&&y<y2)||(y2<=y&&y<y1)){double t=double(y-y1)/double(y2-y1);xs.push_back(static_cast<int>(std::lround(x1+t*(x2-x1))));}}
        std::sort(xs.begin(),xs.end());for(std::size_t i=1;i<xs.size();i+=2)for(int x=xs[i-1];x<=xs[i];++x)im.set(x,y,r,g,b);
    }
}
void write_jpeg(const fs::path&path,const Mesh&m,int cid,int tid,int size){
    Image im(size,size);const int pad=static_cast<int>(0.08*size), top=static_cast<int>(0.12*size);
    const int side=size-2*pad, y0=top+(size-top-pad-side)/2;
    auto pix=[&](P2 p){return std::pair<int,int>{pad+static_cast<int>((p.x+1)*0.5*side),y0+side-static_cast<int>((p.y+1)*0.5*side)};};
    auto draw=[&](const auto&cell,int phase){
        std::vector<std::pair<int,int>> poly;for(int id:cell.v)poly.push_back(pix(m.nodes[id].p));
        if(phase<0)fill_poly(im,poly,210,229,248);else fill_poly(im,poly,245,220,206);
        for(std::size_t i=0;i<poly.size();++i)im.line(poly[i].first,poly[i].second,poly[(i+1)%poly.size()].first,poly[(i+1)%poly.size()].second,65,65,65,1);
    };
    for(const auto&t:m.tris)draw(t,t.phase);for(const auto&q:m.quads)draw(q,q.phase);
    auto[a,b]=std::pair{pix(m.iface.a),pix(m.iface.b)};im.line(a.first,a.second,b.first,b.second,210,30,30,4);
    for(const auto&n:m.nodes)if(n.primary){auto p=pix(n.p);for(int yy=-2;yy<=2;++yy)for(int xx=-2;xx<=2;++xx)if(xx*xx+yy*yy<=4)im.set(p.first+xx,p.second+yy,25,25,25);}
    im.text(pad,18,"CASE "+std::to_string(cid)+" TEMPLATE "+(tid<10?"0":"")+std::to_string(tid),2);
    im.text(pad,size-25,"NEG  POS  TRI  QUAD",1);
    FILE*f=std::fopen(path.string().c_str(),"wb");if(!f)throw std::runtime_error("cannot open JPEG output");
    jpeg_compress_struct c{};jpeg_error_mgr e{};c.err=jpeg_std_error(&e);jpeg_create_compress(&c);jpeg_stdio_dest(&c,f);
    c.image_width=size;c.image_height=size;c.input_components=3;c.in_color_space=JCS_RGB;jpeg_set_defaults(&c);jpeg_set_quality(&c,92,TRUE);jpeg_start_compress(&c,TRUE);
    while(c.next_scanline<c.image_height){JSAMPROW row=&im.px[static_cast<std::size_t>(c.next_scanline*size*3)];jpeg_write_scanlines(&c,&row,1);}
    jpeg_finish_compress(&c);jpeg_destroy_compress(&c);std::fclose(f);
}
void write_text_atomic(const fs::path&final,const std::string&text){
    const auto tmp=final.string()+".tmp";{std::ofstream o(tmp,std::ios::binary);if(!o)throw std::runtime_error("cannot write "+tmp);o<<text;if(!o)throw std::runtime_error("write failed "+tmp);}
    fs::rename(tmp,final);
}
void publish_pair(const fs::path&dir,const std::string&stem,const std::string&dat,const Mesh&m,int cid,int tid,int image_size){
    const fs::path dat_tmp=dir/(stem+".dat.tmp"),jpg_tmp=dir/(stem+".jpg.tmp"),dat_final=dir/(stem+".dat"),jpg_final=dir/(stem+".jpg");
    try{
        {std::ofstream o(dat_tmp,std::ios::binary);if(!o)throw std::runtime_error("cannot open temporary dataset");o<<dat;if(!o)throw std::runtime_error("dataset write failed");}
        write_jpeg(jpg_tmp,m,cid,tid,image_size);
        fs::rename(dat_tmp,dat_final);fs::rename(jpg_tmp,jpg_final);
    }catch(...){
        std::error_code ec;
        fs::remove(dat_tmp,ec);fs::remove(jpg_tmp,ec);
        fs::remove(dat_final,ec);fs::remove(jpg_final,ec);
        throw;
    }
}
std::string manifest_json(const Config&cfg,const GenerationSummary&s,const std::vector<std::tuple<int,int,std::string>>&accepted){
    std::ostringstream o;o<<"{\n  \"schema_version\": \""<<kSchemaVersion<<"\",\n  \"seed\": "<<cfg.seed<<",\n  \"requested\": "<<s.requested<<",\n  \"accepted\": "<<s.accepted<<",\n  \"failed\": "<<s.failed<<",\n  \"templates\": [\n";
    for(std::size_t i=0;i<accepted.size();++i){auto[c,t,p]=accepted[i];o<<"    {\"case\": "<<c<<", \"template\": "<<t<<", \"path\": \""<<p<<"\"}"<<(i+1<accepted.size()?",":"")<<"\n";}
    o<<"  ],\n  \"failures\": [";
    for(std::size_t i=0;i<s.failures.size();++i){std::string x=s.failures[i];for(auto&ch:x)if(ch=='\"')ch='\'';o<<(i?", ":"")<<"\""<<x<<"\"";}
    o<<"]\n}\n";return o.str();
}

} // namespace

int classify_phi(double phi,double tolerance){
    if(!std::isfinite(phi)||!std::isfinite(tolerance)||!(tolerance>0.0))throw std::invalid_argument("phi and tolerance must be finite; tolerance positive");
    if(phi < -tolerance)return -1;if(phi > tolerance)return 1;return 0;
}
void validate_config(const Config&c){
    auto finite=[](double x){return std::isfinite(x);};
    if(c.templates_per_case<=0)throw std::invalid_argument("templates-per-case must be positive");
    if(c.cases.empty())throw std::invalid_argument("cases cannot be empty");
    std::set<int> seen;for(int id:c.cases){if(id!=2&&id!=3&&id!=6&&id!=11)throw std::invalid_argument("supported cases are 2,3,6,11");if(!seen.insert(id).second)throw std::invalid_argument("duplicate case id");}
    if(!finite(c.phi_zero_tol)||!(c.phi_zero_tol>0))throw std::invalid_argument("phi-zero-tol must be finite and positive");
    if(!finite(c.min_edge_fraction)||!(c.min_edge_fraction>0&&c.min_edge_fraction<0.5))throw std::invalid_argument("min-edge-fraction must lie strictly between 0 and 0.5");
    if(!finite(c.target_edge_length)||!(c.target_edge_length>0))throw std::invalid_argument("target-edge-length must be finite and positive");
    if(!finite(c.min_triangle_quality)||c.min_triangle_quality<0||c.min_triangle_quality>1)throw std::invalid_argument("min-triangle-quality must be in [0,1]");
    if(!finite(c.min_quad_quality)||c.min_quad_quality<0||c.min_quad_quality>1)throw std::invalid_argument("min-quad-quality must be in [0,1]");
    if(c.max_nodes==0)throw std::invalid_argument("max-nodes must be positive");
    if(c.max_attempts_per_template<=0)throw std::invalid_argument("max-attempts-per-template must be positive");
    if(c.max_smoothing_passes<0)throw std::invalid_argument("max-smoothing-passes must be nonnegative");
    if(c.image_size<=0)throw std::invalid_argument("image-size must be positive");
}
GenerationSummary generate_dataset(const Config&cfg){
    validate_config(cfg);GenerationSummary s;s.requested=cfg.templates_per_case*static_cast<int>(cfg.cases.size());
    if(fs::exists(cfg.output)&&!cfg.overwrite){
        if(!fs::is_directory(cfg.output))throw std::runtime_error("output path exists and is not a directory");
        const bool nonempty=fs::directory_iterator(cfg.output)!=fs::directory_iterator{};
        if(nonempty)throw std::runtime_error("output directory is not empty; pass --overwrite to replace generated content");
    }
    if(cfg.overwrite&&fs::exists(cfg.output))fs::remove_all(cfg.output);
    fs::create_directories(cfg.output);
    std::vector<std::tuple<int,int,std::string>> accepted;
    for(int cid:cfg.cases){
        const auto dir=cfg.output/("case_"+std::to_string(cid));fs::create_directories(dir);
        std::unordered_set<std::string> signatures;
        for(int ti=0;ti<cfg.templates_per_case;++ti){
            bool ok=false;std::string last;
            for(int attempt=0;attempt<cfg.max_attempts_per_template&&!ok;++attempt){
                try{
                    auto c=build_candidate(cfg,cid,ti,attempt);c.signature=signature_of(c.mesh);
                    if(!signatures.insert(c.signature).second)throw std::runtime_error("duplicate geometry/connectivity dataset");
                    const int out_index=ti+1;std::ostringstream stem;stem<<"temp_"<<std::setw(2)<<std::setfill('0')<<out_index;
                    const auto dat=serialize(c.mesh,cfg,cid,out_index,cfg.seed);
                    validate_serialized_roundtrip(dat,c.mesh);
                    publish_pair(dir,stem.str(),dat,c.mesh,cid,out_index,cfg.image_size);
                    accepted.emplace_back(cid,out_index,(fs::path("case_"+std::to_string(cid))/(stem.str()+".dat")).generic_string());
                    ++s.accepted;ok=true;
                }catch(const std::exception&e){last=e.what();}
            }
            if(!ok){++s.failed;s.failures.push_back("case "+std::to_string(cid)+" template "+std::to_string(ti+1)+": "+last);}
        }
    }
    s.success=(s.accepted==s.requested);
    write_text_atomic(cfg.output/"manifest.json",manifest_json(cfg,s,accepted));
    return s;
}

} // namespace pvmls
