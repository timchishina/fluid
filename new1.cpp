#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "simulator.hpp"

namespace {
#ifdef SIZES
using Simulator = types::Simulator<types::TypesList<TYPES>, SIZES>;
#else
using Simulator = types::Simulator<types::TypesList<TYPES>>;
#endif

using std::string;
using std::string_view;

[[nodiscard]] string get_arg(string_view arg_name,
                             int argc,
                             char **argv,
                             string_view default_value) {
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == arg_name) {
            return argv[i + 1];
        }
    }
    return string(default_value);
}

using Field = std::vector<std::vector<char>>;

[[nodiscard]] size_t parse_size_t(std::string_view str) noexcept {
    size_t n{};
    [[maybe_unused]] bool converted =
        std::from_chars(str.data(), str.data() + str.size(), n).ec == std::errc{};
    assert(converted);
    return n;
}

[[nodiscard]] std::pair<size_t, size_t> read_sizes(std::istream &in) {
    string line_with_sizes;
    std::getline(in, line_with_sizes);
    const size_t sep_pos = line_with_sizes.find(' ');
    const size_t n = parse_size_t(string_view(line_with_sizes).substr(0, sep_pos));
    const size_t m = parse_size_t(string_view(line_with_sizes).substr(sep_pos + 1));
    return {n, m};
}

[[nodiscard]] Field read_field_from_path(const string &path) {
    std::ifstream fin(path);
    fin.exceptions(std::ios::badbit);

    const auto [n, m] = read_sizes(fin);

    Field field{};
    field.reserve(n);

    while (fin) {
        std::string line;
        line.reserve(m + 4);
        std::getline(fin, line);
        if (line.empty()) {
            break;
        }
        line.pop_back();
        assert(line.size() == m);
        field.push_back(std::vector<char>{line.begin(), line.end()});
    }

    assert(field.size() == n);
    return field;
}

}  // namespace

int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();
    string P_type  = get_arg("--p-type", argc, argv, "FAST_FIXED(32,16)");
    string V_type  = get_arg("--v-type", argc, argv, "FIXED(31,17)");
    string Vf_type = get_arg("--v-flow-type", argc, argv, "DOUBLE");
    string path    = get_arg("--path", argc, argv, "./field");

    const Simulator simulator = Simulator::from_params(types::SimulationParams{
        .p_type_name      = P_type,
        .v_type_name      = V_type,
        .v_flow_type_name = Vf_type,
    });
    simulator.start_on_field(types::Context{
        .field = read_field_from_path(path),
    });
     auto end = std::chrono::high_resolution_clock::now();
      // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Output the time taken in seconds
    std::cout << "Time taken: " << duration.count() << " seconds\n";

    return 0;
}
