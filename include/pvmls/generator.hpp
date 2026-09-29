#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace pvmls {

struct Config {
    std::filesystem::path output{"./pvmls_templates"};
    int templates_per_case{15};
    std::vector<int> cases{2, 3, 6, 11};
    std::uint64_t seed{24301ULL};
    double phi_zero_tol{1e-10};
    double min_edge_fraction{1e-3};
    double target_edge_length{0.3};
    double min_triangle_quality{0.05};
    double min_quad_quality{0.2};
    std::size_t max_nodes{512};
    int max_attempts_per_template{100};
    int max_smoothing_passes{50};
    int image_size{1024};
    bool overwrite{false};
};

struct GenerationSummary {
    int requested{0};
    int accepted{0};
    int failed{0};
    bool success{false};
    std::vector<std::string> failures;
};

int classify_phi(double phi, double tolerance);
void validate_config(const Config& config);
GenerationSummary generate_dataset(const Config& config);

} // namespace pvmls
