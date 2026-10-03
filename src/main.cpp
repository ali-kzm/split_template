#include "pvmls/generator.hpp"

#include <charconv>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void usage() {
    std::cout
        << "pvmls-template-generator [options]\n"
        << "  --output <directory>              default ./pvmls_templates\n"
        << "  --templates-per-case <N>          default 15\n"
        << "  --cases <ids>                     comma-separated, default 2,3,6,11\n"
        << "  --seed <integer>                  default 24301\n"
        << "  --phi-zero-tol <real>             default 1e-10\n"
        << "  --min-edge-fraction <real>        default 1e-3\n"
        << "  --target-edge-length <real>       legacy spacing hint, default 0.3\n"
        << "  --target-elements <N>             desired primary cells, default 10\n"
        << "  --max-elements <N>                hard primary-cell cap, default 20\n"
        << "  --triangle-aspect-threshold <r>   keep triangle whole below r, default 3\n"
        << "  --quad-aspect-threshold <r>       regular/thin quad threshold, default 3\n"
        << "  --pentagon-small-edge-fraction <r> tiny-edge trigger, default 0.15\n"
        << "  --min-triangle-quality <real>     default 0.05\n"
        << "  --min-quad-quality <real>         default 0.2\n"
        << "  --max-nodes <N>                   default 512\n"
        << "  --max-attempts-per-template <N>   default 100\n"
        << "  --max-smoothing-passes <N>        default 50\n"
        << "  --image-size <N>                  default 1024\n"
        << "  --overwrite\n";
}

template <class T>
T parse_number(const std::string& text) {
    std::stringstream ss(text);
    T value{};
    ss >> value;
    if (!ss || !ss.eof()) {
        throw std::invalid_argument("invalid numeric value: " + text);
    }
    return value;
}

std::vector<int> parse_cases(const std::string& text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) throw std::invalid_argument("empty case id");
        out.push_back(parse_number<int>(token));
    }
    if (out.empty()) throw std::invalid_argument("--cases cannot be empty");
    return out;
}

} // namespace

int main(int argc, char** argv) {
    try {
        pvmls::Config cfg;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto need = [&](const char* name) -> std::string {
                if (i + 1 >= argc) throw std::invalid_argument(std::string(name) + " requires a value");
                return argv[++i];
            };
            if (arg == "--help" || arg == "-h") { usage(); return 0; }
            else if (arg == "--output") cfg.output = need("--output");
            else if (arg == "--templates-per-case") cfg.templates_per_case = parse_number<int>(need(arg.c_str()));
            else if (arg == "--cases") cfg.cases = parse_cases(need(arg.c_str()));
            else if (arg == "--seed") cfg.seed = parse_number<std::uint64_t>(need(arg.c_str()));
            else if (arg == "--phi-zero-tol") cfg.phi_zero_tol = parse_number<double>(need(arg.c_str()));
            else if (arg == "--min-edge-fraction") cfg.min_edge_fraction = parse_number<double>(need(arg.c_str()));
            else if (arg == "--target-edge-length") cfg.target_edge_length = parse_number<double>(need(arg.c_str()));
            else if (arg == "--target-elements") cfg.target_elements = parse_number<int>(need(arg.c_str()));
            else if (arg == "--max-elements") cfg.max_elements = parse_number<int>(need(arg.c_str()));
            else if (arg == "--triangle-aspect-threshold") cfg.triangle_aspect_threshold = parse_number<double>(need(arg.c_str()));
            else if (arg == "--quad-aspect-threshold") cfg.quad_aspect_threshold = parse_number<double>(need(arg.c_str()));
            else if (arg == "--pentagon-small-edge-fraction") cfg.pentagon_small_edge_fraction = parse_number<double>(need(arg.c_str()));
            else if (arg == "--min-triangle-quality") cfg.min_triangle_quality = parse_number<double>(need(arg.c_str()));
            else if (arg == "--min-quad-quality") cfg.min_quad_quality = parse_number<double>(need(arg.c_str()));
            else if (arg == "--max-nodes") cfg.max_nodes = parse_number<std::size_t>(need(arg.c_str()));
            else if (arg == "--max-attempts-per-template") cfg.max_attempts_per_template = parse_number<int>(need(arg.c_str()));
            else if (arg == "--max-smoothing-passes") cfg.max_smoothing_passes = parse_number<int>(need(arg.c_str()));
            else if (arg == "--image-size") cfg.image_size = parse_number<int>(need(arg.c_str()));
            else if (arg == "--overwrite") cfg.overwrite = true;
            else throw std::invalid_argument("unknown argument: " + arg);
        }

        pvmls::validate_config(cfg);
        const auto summary = pvmls::generate_dataset(cfg);
        std::cout << "accepted " << summary.accepted << "/" << summary.requested
                  << " templates; failed=" << summary.failed << "\n";
        for (const auto& f : summary.failures) std::cerr << f << "\n";
        return summary.success ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
