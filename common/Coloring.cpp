#include <algorithm>
#include <limits>
#include <numeric>
#include "Types.hpp"
#include "Coloring.hpp"

// グラフ G に対して貪欲法で彩色を行う
// 戻り値: 使用した色の数
// 色情報を格納した配列 color も返す
// 彩色方針は「次数降順＋ID昇順タイブレーク」
int Greedy_Coloring(const Graph& G, std::vector<int>& color)
{
    const std::size_t N = G.n;
    color.assign(N, -1); // -1: uncolored

    struct VertexDegree {
        int v;
        std::size_t deg;
    };

    std::vector<VertexDegree> order;
    order.reserve(N);
    for (int v = 0; v < static_cast<int>(N); ++v) {
        order.push_back({v, G.adj[v].size()});
    }

    // 次数降順、同次数ならID昇順
    std::sort(order.begin(), order.end(),
              [](const VertexDegree& a, const VertexDegree& b) {
                  if (a.deg != b.deg) return a.deg > b.deg; // 降順
                  return a.v < b.v;                         // タイブレーク：ID昇順
              });

    // 近傍色のマーキング（タイムスタンプ法）
    std::vector<int> mark(N + 1, -1); // 色番号は最大でも N-1
    int stamp = 0;
    int max_color = -1;

    for (const auto& vd : order) {
        int u = vd.v;
        ++stamp;

        for (const auto& e : G.adj[u]) {
            int c = color[e.to];
            if (c >= 0) mark[c] = stamp;
        }

        // 最小許容色を割り当て
        int c = 0;
        while (c <= max_color && mark[c] == stamp) ++c;
        if (c == max_color + 1) ++max_color;
        color[u] = c;
    }

    return max_color + 1; // 使用色数
}

// グラフ G に対して貪欲法で彩色を行う
// ただし、使用可能な既存色の中で、色クラスサイズが最小の色を優先する
// 新しい色は、使用可能な既存色が存在しない場合のみ作る
// 戻り値: 使用した色の数
int Greedy_Coloring_Balanced(const Graph& G, std::vector<int>& color)
{
    const std::size_t N = G.n;
    color.assign(N, -1); // -1: uncolored
    if (N == 0) return 0;

    struct VertexDegree {
        int v;
        std::size_t deg;
    };

    std::vector<VertexDegree> order;
    order.reserve(N);
    for (int v = 0; v < static_cast<int>(N); ++v) {
        order.push_back({v, G.adj[v].size()});
    }

    // 次数降順、同次数ならID昇順
    std::sort(order.begin(), order.end(),
              [](const VertexDegree& a, const VertexDegree& b) {
                  if (a.deg != b.deg) return a.deg > b.deg;
                  return a.v < b.v;
              });

    // 近傍色のマーキング（タイムスタンプ法）
    std::vector<int> mark(N + 1, -1);
    std::vector<int> color_size(N + 1, 0);
    int stamp = 0;
    int max_color = -1;

    for (const auto& vd : order) {
        const int u = vd.v;
        ++stamp;

        for (const auto& e : G.adj[u]) {
            const int c = color[e.to];
            if (c >= 0) mark[c] = stamp;
        }

        int best_color = -1;
        int best_size = std::numeric_limits<int>::max();
        for (int c = 0; c <= max_color; ++c) {
            if (mark[c] != stamp && color_size[c] < best_size) {
                best_color = c;
                best_size = color_size[c];
            }
        }

        // 使用可能な既存色がない場合のみ、新しい色を作る
        if (best_color == -1) {
            best_color = ++max_color;
        }

        color[u] = best_color;
        ++color_size[best_color];
    }

    return max_color + 1;
}

// Greedy_Coloring で色数 nc を確定した後、nc 色を固定したまま色クラスサイズを均等化する
// 新しい色は作らず、隣接頂点と同じ色にならない合法な頂点移動のみを行う
// 戻り値: 使用色数 nc
int Greedy_Coloring_FixedNc_Balanced(const Graph& G, std::vector<int>& color)
{
    const int N = static_cast<int>(G.n);
    const int nc = Greedy_Coloring(G, color);
    if (nc <= 1) return nc;

    std::vector<int> color_size(nc, 0);
    for (int v = 0; v < N; ++v) {
        ++color_size[color[v]];
    }

    auto can_move_to_color = [&](int u, int c_to) {
        for (const auto& e : G.adj[u]) {
            if (color[e.to] == c_to) return false;
        }
        return true;
    };

    std::vector<int> from_order(nc);
    std::vector<int> to_order(nc);
    const int max_passes = 10 * N;

    for (int pass = 0; pass < max_passes; ++pass) {
        std::iota(from_order.begin(), from_order.end(), 0);
        std::iota(to_order.begin(), to_order.end(), 0);

        std::sort(from_order.begin(), from_order.end(),
                  [&](int a, int b) {
                      if (color_size[a] != color_size[b]) {
                          return color_size[a] > color_size[b];
                      }
                      return a < b;
                  });
        std::sort(to_order.begin(), to_order.end(),
                  [&](int a, int b) {
                      if (color_size[a] != color_size[b]) {
                          return color_size[a] < color_size[b];
                      }
                      return a < b;
                  });

        bool moved = false;
        for (int c_from : from_order) {
            for (int c_to : to_order) {
                if (c_from == c_to ||
                    color_size[c_from] <= color_size[c_to] + 1) {
                    continue;
                }

                for (int u = 0; u < N; ++u) {
                    if (color[u] != c_from || !can_move_to_color(u, c_to)) {
                        continue;
                    }

                    color[u] = c_to;
                    --color_size[c_from];
                    ++color_size[c_to];
                    moved = true;
                    break;
                }

                if (moved) break;
            }
            if (moved) break;
        }

        if (!moved) break;
    }

    return nc;
}

// 頻度順に色ラベルを付け替える
// 同数なら旧ラベルの昇順を優先（安定なタイブレーク）
void RelabelColorsByClassSize(std::vector<int>& color) {
    if (color.empty()) return;

    int nc = 1 + *std::max_element(color.begin(), color.end());
    std::vector<int> cnt(nc, 0);
    for (int c : color) if (c >= 0) ++cnt[c];

    std::vector<int> order(nc);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){
                  if (cnt[a] != cnt[b]) return cnt[a] > cnt[b]; // 大きい順
                  return a < b;                                  // タイブレーク
              });

    // old -> new の写像を作る
    std::vector<int> new_id(nc, -1);
    for (int i = 0; i < nc; ++i) new_id[ order[i] ] = i;

    // 配列を書き換え
    for (int& c : color) if (c >= 0) c = new_id[c];
}
