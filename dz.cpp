#include "fixed.hpp"
#include <array>
#include <map>
#include <tuple>
#include <random>
#include <cassert>
#include <limits>
#include <algorithm> // For std::upper_bound

// Existing code in dz1.cpp
template<typename T, size_t N, size_t M>
class SimulationFramework {
public:
    using FixedType = Fixed<32, 16>; // Example usage of Fixed<N, K>

    class VectorField {
        std::map<std::tuple<int, int, int, int>, T> field;
    public:
        T get(int x, int y, int dx, int dy) {
            return field[{x, y, dx, dy}];
        }

        void add(int x, int y, int dx, int dy, T value) {
            field[{x, y, dx, dy}] += value;
        }
    };

    class ParticleParams {
    public:
        void swap_with(int x, int y) {
            // Implement the logic for swapping particles at position (x, y)
        }
    };

    class Simulation {
    public:
        Simulation() {
            initializeField();
            initializeParameters();
        }

        void run() {
            for (size_t i = 0; i < maxTicks; ++i) {
                update();
                if (propagate()) {
                    printField(i);
                }
            }
        }

    private:
        static constexpr size_t maxTicks = 1'000'000;
        static constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

        char field[N][M + 1];
        T p[N][M]{}, old_p[N][M];
        VectorField velocity{}, velocity_flow{};
        int last_use[N][M]{};
        int UT = 0;
        std::mt19937 rnd{1337};
        int dirs[N][M]{};

        void initializeField() {
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    field[i][j] = (i == 0 || i == N-1 || j == 0 || j == M-1) ? '#' : '.';
                }
                field[i][M] = '\0';
            }
        }

        void initializeParameters() {
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    p[i][j] = T(0);
                    old_p[i][j] = T(0);
                    last_use[i][j] = 0;
                    dirs[i][j] = -1;
                }
            }
        }

        void update() {
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    old_p[i][j] = p[i][j];
                    p[i][j] += T(1); // Example update logic
                }
            }
        }

        bool propagate() {
            bool changed = false;
            for (size_t i = 1; i < N-1; ++i) {
                for (size_t j = 1; j < M-1; ++j) {
                    if (field[i][j] == '.') {
                        T new_value = (p[i-1][j] + p[i+1][j] + p[i][j-1] + p[i][j+1]) / T(4);
                        if (new_value != p[i][j]) {
                            p[i][j] = new_value;
                            changed = true;
                        }
                    }
                }
            }
            return changed;
        }

        void printField(size_t tick) {
            std::cout << "Tick " << tick << ":\n";
            for (size_t x = 0; x < N; ++x) {
                std::cout << field[x] << "\n";
            }
        }

        std::tuple<T, bool, std::pair<int, int>> propagate_flow(int x, int y, T lim) {
            last_use[x][y] = UT - 1;
            T ret = 0;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                    auto cap = velocity.get(x, y, dx, dy);
                    auto flow = velocity_flow.get(x, y, dx, dy);
                    if (flow == cap) {
                        continue;
                    }
                    auto vp = std::min(lim, cap - flow);
                    if (last_use[nx][ny] == UT - 1) {
                        velocity_flow.add(x, y, dx, dy, vp);
                        last_use[x][y] = UT;
                        return {vp, 1, {nx, ny}};
                    }
                    auto [t, prop, end] = propagate_flow(nx, ny, vp);
                    ret += t;
                    if (prop) {
                        velocity_flow.add(x, y, dx, dy, t);
                        last_use[x][y] = UT;
                        return {t, prop && end != std::pair(x, y), end};
                    }
                }
            }
            last_use[x][y] = UT;
            return {ret, 0, {0, 0}};
        }

        T random01() {
            return T::from_raw((rnd() & ((1 << 16) - 1)));
        }

        void propagate_stop(int x, int y, bool force = false) {
            if (!force) {
                bool stop = true;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                        stop = false;
                        break;
                    }
                }
                if (!stop) {
                    return;
                }
            }
            last_use[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                    continue;
                }
                propagate_stop(nx, ny);
            }
        }

        T move_prob(int x, int y) {
            T sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    continue;
                }
                sum += v;
            }
            return sum;
        }

        T rho[256];
        bool propagate_move(int x, int y, bool is_first) {
            last_use[x][y] = UT - is_first;
            bool ret = false;
            int nx = -1, ny = -1;
            do {
                std::array<T, deltas.size()> tres;
                T sum = 0;
                for (size_t i = 0; i < deltas.size(); ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                        tres[i] = sum;
                        continue;
                    }
                    auto v = velocity.get(x, y, dx, dy);
                    if (v < 0) {
                        tres[i] = sum;
                        continue;
                    }
                    sum += v;
                    tres[i] = sum;
                }

                if (sum == 0) {
                    break;
                }

                T p = random01() * sum;
                size_t d = std::upper_bound(tres.begin(), tres.end(), p) - tres.begin();

                auto [dx, dy] = deltas[d];
                nx = x + dx;
                ny = y + dy;
                assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

                ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
            } while (!ret);
            last_use[x][y] = UT;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                    propagate_stop(nx, ny);
                }
            }
            if (ret) {
                if (!is_first) {
                    ParticleParams pp{};
                    pp.swap_with(x, y);
                    pp.swap_with(nx, ny);
                    pp.swap_with(x, y);
                }
            }
            return ret;
        }
    };
};

int main() {
    SimulationFramework<int, 36, 84>::Simulation sim;
    sim.run();
    return 0;
}