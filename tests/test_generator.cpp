#include "pvmls/generator.hpp"

extern "C" {
#include <jpeglib.h>
}

#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static std::string slurp(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

static void assert_jpeg(const fs::path& p, int expected_size) {
    FILE* f = std::fopen(p.string().c_str(), "rb");
    assert(f);
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
    assert(s.success);
    assert(s.accepted == 4);
    for (int id : cfg.cases) {
        const auto dat = cfg.output / ("case_" + std::to_string(id)) / "temp_01.dat";
        const auto jpg = cfg.output / ("case_" + std::to_string(id)) / "temp_01.jpg";
        assert(fs::exists(dat));
        assert(fs::exists(jpg));
        const auto text = slurp(dat);
        assert(text.find("VALIDATION OK") != std::string::npos);
        assert(text.find("PVMLS_TEMPLATE_DATA 1.0.0") != std::string::npos);
        assert_jpeg(jpg, cfg.image_size);
    }

    // Case 11 must keep fixed interface geometry while varying interior layout.
    cfg.output = root / "case11_variation";
    cfg.templates_per_case = 2;
    cfg.cases = {11};
    cfg.seed = 2002;
    cfg.overwrite = true;
    s = pvmls::generate_dataset(cfg);
    assert(s.success && s.accepted == 2);
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
    assert(pvmls::generate_dataset(r1).success);
    assert(pvmls::generate_dataset(r2).success);
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
    assert(near_s.success);

    fs::remove_all(root);
    std::cout << "all tests passed\n";
    return 0;
}
