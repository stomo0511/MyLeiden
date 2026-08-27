#include <algorithm>
#include <unordered_map>

#include "Types.hpp"
#include "mm_io.hpp"
#include "BlockIO.hpp"

// ---- 出力：ブロック情報（1始まり）
// 1行目: <#blocks>
// 2行目以降: "<node_id> <block_id>" (node_id, block_id は 1始まり)
void WriteBlockInfo_1Based(const std::vector<int>& block_of, const std::string& out_path) {
    std::ofstream fout(out_path);
    if (!fout)
        throw std::runtime_error("cannot open output: " + out_path);

    // ブロック数を算出
    int nb = *std::max_element(block_of.begin(), block_of.end()) + 1;

    // 先頭行にブロック数を出力
    fout << nb << "\n";

    for (int i = 0; i < (int)block_of.size(); ++i) {
        int bid = block_of[i];
        if (bid < 0)
            throw std::runtime_error("unassigned node in block_of");
        fout << (i+1) << " " << (bid+1) << "\n";
    }
    // for (int old = 0; old < part.n; ++old) {
    //     int bid = part.block_of[old];
    //     if (bid < 0) throw std::runtime_error("unassigned node in block_of");
    //     fout << (old+1) << " " << (bid+1) << "\n";
    // }
}

// ---- 出力：ブロック色情報（1始まり）
// 1行目: <#colors>
// 2行目以降: "<block_id> <color_id>" (block_id, color_id は 1始まり)
void WriteBlockColor_1Based(const std::vector<int>& block_color, int nc, const std::string& out_path) {
    std::ofstream fout(out_path);
    if (!fout)
        throw std::runtime_error("cannot open output: " + out_path);

    fout << nc << "\n";
    for (int bid = 0; bid < (int)block_color.size(); ++bid) {
        int c = block_color[bid];
        if (!(0 <= c && c < nc)) throw std::runtime_error("block_color out of range");
        fout << (bid+1) << " " << (c+1) << "\n";
    }
}

// ブロックグラフの作成
Graph BuildBlockGraph(const Graph& G,
                      const std::vector<int>& block_of,
                      BlockEdgeWeight mode = BlockEdgeWeight::Binary)
{
    const int nb = *std::max_element(block_of.begin(), block_of.end()) + 1;
    Graph T(nb);

    // 集計用: (min(bk,bl), max(bk,bl)) -> weight
    struct PairHash {
        size_t operator()(const std::pair<int,int>& p) const {
            return (static_cast<size_t>(p.first) << 32) ^ static_cast<size_t>(p.second);
        }
    };
    std::unordered_map<std::pair<int,int>, double, PairHash> agg;
    agg.reserve(num_edges(G));

    for_each_undirected_edge(G, [&](int u, int v, double w) {
        int bk = block_of[u];
        int bl = block_of[v];
        if (bk == bl) return;
        if (bk > bl) std::swap(bk, bl);

        double inc = 0.0;
        switch (mode) {
            case BlockEdgeWeight::Binary: inc = 1.0; break;
            case BlockEdgeWeight::Count:  inc = 1.0; break;
            case BlockEdgeWeight::AbsSum: inc = std::abs(w); break;
        }
        auto key = std::make_pair(bk, bl);
        if (mode == BlockEdgeWeight::Binary) {
            // 1 回でも見つかれば 1 にして固定
            auto it = agg.find(key);
            if (it == agg.end()) agg.emplace(key, 1.0);
        } else {
            agg[key] += inc;
        }
    });

    // 集計結果を T に反映
    for (const auto& kv : agg) {
        int bk = kv.first.first;
        int bl = kv.first.second;
        add_undirected_edge(T, bk, bl, (mode == BlockEdgeWeight::Binary) ? 1.0 : kv.second);
    }
    return T;
}

/**
 * .blk ファイル（1-based）を読み込み、0-based の vector に変換する
 */
int ReadBlockInfo_1Based(const std::string& path, int N, std::vector<int>& block_of) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open block file: " + path);

    int nb;
    fin >> nb; // 1行目: ブロック数
    block_of.assign(N, -1);

    int u_raw, b_raw;
    while (fin >> u_raw >> b_raw) {
        int u = u_raw - 1; // 1-based -> 0-based
        int b = b_raw - 1;
        if (u >= 0 && u < N) block_of[u] = b;
    }
    return nb;
}

/**
 * .bcol ファイル（1-based）を読み込み、0-based の vector に変換する
 */
int ReadBlockColor_1Based(const std::string& path, int nb, std::vector<int>& block_color) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open color file: " + path);

    int nc;
    fin >> nc; // 1行目: 彩色数
    block_color.assign(nb, -1);

    int b_raw, c_raw;
    while (fin >> b_raw >> c_raw) {
        int b = b_raw - 1; // 1-based -> 0-based
        int c = c_raw - 1;
        if (b >= 0 && b < nb) block_color[b] = c;
    }
    return nc;
}
