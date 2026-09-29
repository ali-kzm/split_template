#include "pvmls/generator.hpp"

extern "C" {
#include <jpeglib.h>
}

#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

static std::string slurp(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

static std::tuple<double,double,double,double,double,double> read_interface(const fs::path& p) {
    std::ifstream in(p);
    std::string line;
    double x0=0.0,y0=0.0,t0=0.0,x1=0.0,y1=0.0,t1=0.0;
    bool have0=false,have1=false;
    while (std::getline(in,line)) {
        if (line.rfind("INTERFACE_ENDPOINT ",0)!=0) continue;
        std::istringstream row(line);
        std::string tag, phi_tag, edge_tag, parameter_tag;
        int endpoint=-1, edge=-1;
        double x=0.0,y=0.0,phi=1.0,t=0.0;
        row >> tag >> endpoint >> x >> y >> phi_tag >> phi >> edge_tag >> edge >> parameter_tag >> t;
        assert(row);
        assert(phi_tag=="phi");
        assert(phi==0.0);
        assert(edge_tag=="edge");
        assert(parameter_tag=="parameter");
        if (endpoint==0) { x0=x; y0=y; t0=t; have0=true; }
        if (endpoint==1) { x1=x; y1=y; t1=t; have1=true; }
    }
    assert(have0&&have1);
    return {x0,y0,t0,x1,y1,t1};
}

static int stratum_of(double t, int n, double lo, double hi) {
    const double u=(t-lo)/(hi-lo);
    return std::clamp(static_cast<int>(std::floor(u*n)),0,n-1);
}

static int read_section_count(const fs::path& p, const std::string& prefix) {
    std::ifstream in(p);
    std::string line;
    while (std::getline(in,line)) {
        if (line.rfind(prefix,0)!=0) continue;
        std::istringstream row(line);
        std::string tag;
        int n=-1;
        row >> tag >> n;
        assert(row && n>=0);
        return n;
    }
    assert(false);
    return -1;
}

static int read_element_count(const fs::path& p) {
    std::ifstream in(p);
    std::string line;
    while (std::getline(in,line)) {
        if (line.rfind("METRICS ",0)!=0) continue;
        const auto pos=line.find(" element_count ");
        assert(pos!=std::string::npos);
        std::istringstream row(line.substr(pos + std::string(" element_count ").size()));
        int n=-1;
        row >> n;
        assert(row && n>=0);
        return n;
    }
    assert(false);
    return -1;
}

static void require_success(const pvmls::GenerationSummary& s) {
    if (s.success) return;
    std::cerr << "generation failed: accepted=" << s.accepted << "/" << s.requested << "\n";
    for (const auto& f : s.failures) std::cerr << "  " << f << "\n";
    std::abort();
}

static void assert_jpeg(const fs::path& p, int expected_size) {
    FILE* f = std::fopen(p.string().c_str(), "rb");
    if (!f) {
        std::cerr << "missing JPEG: " << p << "\n";
        std::abort();
    }
    jpeg_decompress_struct cinfo{};
    jpeg_error_mgr jerr{};
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, f);
    assert(jpeg_read_header(&cinfo, TRUE) == JPEG_HEADER_OK);
    assert(static_cast<int>(cinfo.image_width) == expected_size);
    assert(static_cast<int>(cinfo.image_height) == expected_size);
    jpeg_destroy_decompress(&cinfo);
    std::fclose(f);
}

int main() {
    assert(pvmls::classify_phi(-1.1e-10, 1e-10) == -1);
    assert(pvmls::classify_phi(-1.0e-10, 1e-10) == 0);
    assert(pvmls::classify_phi( 1.0e-10, 1e-10) == 0);
    assert(pvmls::classify_phi( 1.1e-10, 1e-10) == 1);

    const fs::path root = fs::temp_directory_path() / "pvmls_generator_tests";
    fs::remove_all(root);

    pvmls::Config cfg;
    cfg.output = root / "families";
    cfg.templates_per_case = 1;
    cfg.cases = {2, 3, 6, 11};
    cfg.seed = 1001;
    cfg.image_size = 192;
    cfg.overwrite = true;
    auto s = pvmls::generate_dataset(cfg);
    require_success(s);
    assert(s.accepted == 4);
    for (int id : cfg.cases) {
        const auto dat = cfg.output / ("case_" + std::to_string(id)) / "temp_01.dat";
        const auto jpg = cfg.output / ("case_" + std::to_string(id)) / "temp_01.jpg";
        assert(fs::exists(dat));
        assert(fs::exists(jpg));
        const auto text = slurp(dat);
        assert(text.find("VALIDATION OK") != std::string::npos);
        assert(text.find("PVMLS_TEMPLATE_DATA 2.0.0") != std::string::npos);
        const int neg_edges=read_section_count(dat,"INTERFACE_EDGES_NEGATIVE ");
        const int pos_edges=read_section_count(dat,"INTERFACE_EDGES_POSITIVE ");
        assert(neg_edges>=1 && pos_edges>=1);
        assert(neg_edges!=pos_edges);
        assert(read_element_count(dat)<=cfg.max_elements);
        assert_jpeg(jpg, cfg.image_size);
    }

    // Sampling regression: case 2 must contain multiple interface slopes,
    // and case 3 must not contain vertical-flip-equivalent stratum pairs.
    pvmls::Config sampling;
    sampling.output = root / "sampling_regression";
    sampling.templates_per_case = 9;
    sampling.cases = {2, 3};
    sampling.seed = 24301;
    sampling.min_edge_fraction = 0.02;
    sampling.target_edge_length = 0.5;
    sampling.image_size = 96;
    sampling.overwrite = true;
    auto sampling_s = pvmls::generate_dataset(sampling);
    require_success(sampling_s);
    assert(sampling_s.accepted == 18);

    std::set<long long> case2_angles;
    for (int i=1;i<=sampling.templates_per_case;++i) {
        std::ostringstream name; name << "temp_" << std::setw(2) << std::setfill('0') << i << ".dat";
        const auto dat=sampling.output / "case_2" / name.str();
        auto [x0,y0,t0,x1,y1,t1] = read_interface(dat);
        assert(read_section_count(dat,"INTERFACE_EDGES_NEGATIVE ") !=
               read_section_count(dat,"INTERFACE_EDGES_POSITIVE "));
        assert(read_element_count(dat)<=sampling.max_elements);
        (void)t0; (void)t1;
        const double angle=std::atan2(y1-y0,x1-x0);
        case2_angles.insert(static_cast<long long>(std::llround(angle*1e6)));
    }
    assert(case2_angles.size() >= 4);

    std::vector<std::pair<int,int>> case3_strata;
    for (int i=1;i<=sampling.templates_per_case;++i) {
        std::ostringstream name; name << "temp_" << std::setw(2) << std::setfill('0') << i << ".dat";
        auto [x0,y0,t0,x1,y1,t1] = read_interface(sampling.output / "case_3" / name.str());
        (void)x0; (void)y0; (void)x1; (void)y1;
        case3_strata.push_back({
            stratum_of(t0,sampling.templates_per_case,sampling.min_edge_fraction,1.0-sampling.min_edge_fraction),
            stratum_of(t1,sampling.templates_per_case,sampling.min_edge_fraction,1.0-sampling.min_edge_fraction)
        });
    }
    for (std::size_t i=0;i<case3_strata.size();++i) {
        for (std::size_t j=i+1;j<case3_strata.size();++j) {
            const bool vertical_flip =
                case3_strata[i].first + case3_strata[j].first == sampling.templates_per_case-1 &&
                case3_strata[i].second + case3_strata[j].second == sampling.templates_per_case-1;
            assert(!vertical_flip);
        }
    }

    // Case 11 must keep fixed interface geometry while varying interior layout.
    cfg.output = root / "case11_variation";
    cfg.templates_per_case = 2;
    cfg.cases = {11};
    cfg.seed = 2002;
    cfg.overwrite = true;
    s = pvmls::generate_dataset(cfg);
    require_success(s);
    assert(s.accepted == 2);
    const auto a = slurp(cfg.output / "case_11/temp_01.dat");
    const auto b = slurp(cfg.output / "case_11/temp_02.dat");
    assert(a != b);
    assert(a.find("INTERFACE_ENDPOINT 0 -1 -1") != std::string::npos);
    assert(a.find("INTERFACE_ENDPOINT 1 1 1") != std::string::npos);
    assert(b.find("INTERFACE_ENDPOINT 0 -1 -1") != std::string::npos);
    assert(b.find("INTERFACE_ENDPOINT 1 1 1") != std::string::npos);

    // Reproducibility: same seed/config produces byte-identical data files.
    pvmls::Config r1 = cfg;
    r1.templates_per_case = 1;
    r1.cases = {3};
    r1.output = root / "repro_a";
    r1.overwrite = true;
    auto r2 = r1;
    r2.output = root / "repro_b";
    const auto rs1=pvmls::generate_dataset(r1); require_success(rs1);
    const auto rs2=pvmls::generate_dataset(r2); require_success(rs2);
    assert(slurp(r1.output / "case_3/temp_01.dat") == slurp(r2.output / "case_3/temp_01.dat"));

    // Resource-limit failure must be explicit and must not publish an invalid pair.
    pvmls::Config tiny;
    tiny.output = root / "resource_failure";
    tiny.templates_per_case = 1;
    tiny.cases = {2};
    tiny.max_nodes = 8;
    tiny.max_attempts_per_template = 2;
    tiny.image_size = 128;
    tiny.overwrite = true;
    auto fail = pvmls::generate_dataset(tiny);
    assert(!fail.success);
    assert(fail.failed == 1);
    assert(!fs::exists(tiny.output / "case_2/temp_01.dat"));

    // Near-corner cuts are exercised by the first stratum at the configured limit.
    pvmls::Config near;
    near.output = root / "near_corner";
    near.templates_per_case = 2;
    near.cases = {2, 6};
    near.min_edge_fraction = 1e-3;
    near.image_size = 128;
    near.overwrite = true;
    auto near_s = pvmls::generate_dataset(near);
    require_success(near_s);

    fs::remove_all(root);
    std::cout << "all tests passed\n";
    return 0;
}
