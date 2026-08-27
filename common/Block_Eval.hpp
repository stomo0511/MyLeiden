#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "Types.hpp"
#include "Coloring.hpp"
#include "BlockIO.hpp"

/**
 * Block_Eval.hpp
 *
 * Utility routines for evaluating a block partition generated from an
 * adjacency graph.  The routines are kept header-only by using inline
 * functions, so that programs such as eblock.cpp and future iterative
 * Leiden/CPM blocking programs can share the same evaluation code.
 *
 * Assumptions:
 *   - block_of[v] is a 0-based block index.
 *   - block_color[b] is a 0-based color index.
 *   - nb is the number of blocks.
 *   - nc is the number of colors.
 */

struct BlockEvaluationResult {
    int num_vertices = 0;
    int num_edges = 0;
    int num_blocks = 0;
    int num_colors = 0;

    std::vector<int> nodes_per_block;
    std::vector<int> internal_edges_per_block;
    std::vector<int> block_degrees;
    std::vector<int> blocks_per_color;
    std::vector<int> nodes_per_color;

    int internal_edges = 0;
    int external_edges = 0;
    double internal_ratio = 0.0;

    double avg_nodes_per_block = 0.0;
    int max_nodes_per_block = 0;
    int min_nodes_per_block = 0;
    double block_size_stddev = 0.0;
    double block_size_gini = 0.0;

    double mean_internal_avg_degree = 0.0;
    double avg_block_density = 0.0;
    int block_graph_edges = 0;
    double block_graph_density = 0.0;
    double avg_block_graph_degree = 0.0;
    double locality_score = 0.0;

    double modularity_unweighted = 0.0;
    double modularity_weighted = 0.0;
};

inline bool IsValidBlockIndex(int b, int nb)
{
    return 0 <= b && b < nb;
}

inline bool IsValidColorIndex(int c, int nc)
{
    return 0 <= c && c < nc;
}

// ノード数
inline std::vector<int>
CountNodesPerBlock(const std::vector<int>& block_of, int nb)
{
    std::vector<int> cnt(nb, 0);

    for(int b : block_of){
        if(IsValidBlockIndex(b, nb))
            cnt[b] += 1;
    }

    return cnt;
}

// 内部エッジ数
inline std::vector<int>
CountInternalEdges(const Graph& G,
                   const std::vector<int>& block_of,
                   int nb)
{
    std::vector<int> internal(nb, 0);

    for_each_undirected_edge(G, [&](int u, int v, double) {
        int bu = block_of[u];
        int bv = block_of[v];

        if(IsValidBlockIndex(bu, nb) && bu == bv)
            internal[bu] += 1;
    });

    return internal;
}

// ブロックグラフ上の次数 (Binary モデル)
inline std::vector<int>
BlockDegreesBinary(const Graph& T)
{
    int nb = static_cast<int>(num_vertices(T));
    std::vector<int> deg(nb, 0);

    for(int b = 0; b < nb; ++b)
        deg[b] = static_cast<int>(degree(b, T));

    return deg;
}

/**
 * Gini係数を計算する。
 *
 * values は各ブロックのノード数を想定する。
 * 0に近いほど均一、1に近いほど不均一である。
 */
inline double GiniCoefficient(const std::vector<int>& values)
{
    if(values.empty())
        return 0.0;

    std::vector<int> sorted = values;
    std::sort(sorted.begin(), sorted.end());

    double sum = std::accumulate(sorted.begin(), sorted.end(), 0.0);

    if(sum == 0.0)
        return 0.0;

    const int n = static_cast<int>(sorted.size());
    double weighted_sum = 0.0;

    for(int i = 0; i < n; ++i)
        weighted_sum += static_cast<double>(i + 1) * sorted[i];

    return (2.0 * weighted_sum) / (static_cast<double>(n) * sum)
           - static_cast<double>(n + 1) / static_cast<double>(n);
}

inline double Mean(const std::vector<int>& values)
{
    if(values.empty())
        return 0.0;

    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / static_cast<double>(values.size());
}

inline double StddevPopulation(const std::vector<int>& values)
{
    if(values.empty())
        return 0.0;

    double mean = Mean(values);
    double var = 0.0;

    for(int x : values)
        var += (x - mean) * (x - mean);

    var /= static_cast<double>(values.size());

    return std::sqrt(var);
}

inline void CountEdgeLocality(const Graph& G,
                              const std::vector<int>& block_of,
                              int nb,
                              int& internal,
                              int& external)
{
    internal = 0;
    external = 0;

    for_each_undirected_edge(G, [&](int u, int v, double) {
        int bu = block_of[u];
        int bv = block_of[v];

        if(IsValidBlockIndex(bu, nb) && bu == bv)
            internal += 1;
        else
            external += 1;
    });
}

inline double InternalEdgeRatio(int internal, int external)
{
    int total = internal + external;

    if(total == 0)
        return 0.0;

    return static_cast<double>(internal) / static_cast<double>(total);
}

inline double LocalityScore(int internal, int external)
{
    return static_cast<double>(internal) / static_cast<double>(external + 1);
}

inline double AverageBlockDensity(const std::vector<int>& nodes_per_block,
                                  const std::vector<int>& internal_edges)
{
    if(nodes_per_block.empty())
        return 0.0;

    double sum_density = 0.0;

    for(std::size_t b = 0; b < nodes_per_block.size(); ++b){
        int n = nodes_per_block[b];

        if(n <= 1)
            continue;

        double max_edges = static_cast<double>(n) * (n - 1) / 2.0;
        double density = static_cast<double>(internal_edges[b]) / max_edges;

        sum_density += density;
    }

    // 空ブロックや1頂点ブロックも含めて、ブロック数で割る。
    // 旧 eblock.cpp の定義を維持する。
    return sum_density / static_cast<double>(nodes_per_block.size());
}

inline double MeanInternalAverageDegree(const std::vector<int>& nodes_per_block,
                                        const std::vector<int>& internal_edges)
{
    if(nodes_per_block.empty())
        return 0.0;

    double sum = 0.0;

    for(std::size_t b = 0; b < nodes_per_block.size(); ++b){
        if(nodes_per_block[b] > 0){
            double local_avg_deg =
                (2.0 * internal_edges[b]) / nodes_per_block[b];
            sum += local_avg_deg;
        }
    }

    // 空ブロックも含めて、ブロック数で割る。
    // 旧 eblock.cpp の定義を維持する。
    return sum / static_cast<double>(nodes_per_block.size());
}

inline double BlockGraphDensity(const Graph& T)
{
    int nb = static_cast<int>(num_vertices(T));

    if(nb <= 1)
        return 0.0;

    double max_edges = static_cast<double>(nb) * (nb - 1) / 2.0;
    return static_cast<double>(num_edges(T)) / max_edges;
}

// モジュラリティ (未加重)
inline double Modularity_Unweighted(const Graph& G,
                                    const std::vector<int>& block_of)
{
    if(block_of.empty())
        return 0.0;

    int nb = 1 + *std::max_element(block_of.begin(), block_of.end());
    if(nb <= 0)
        return 0.0;

    const double m = static_cast<double>(num_edges(G));
    if(m == 0.0)
        return 0.0;

    std::vector<int> Kb(nb, 0);
    std::vector<int> Lb(nb, 0);

    for(int u = 0; u < static_cast<int>(num_vertices(G)); ++u){
        int b = block_of[u];

        if(IsValidBlockIndex(b, nb))
            Kb[b] += static_cast<int>(degree(u, G));
    }

    for_each_undirected_edge(G, [&](int u, int v, double) {
        int bu = block_of[u];
        int bv = block_of[v];

        if(IsValidBlockIndex(bu, nb) && bu == bv)
            Lb[bu] += 1;
    });

    double Q = 0.0;

    for(int b = 0; b < nb; ++b){
        double L = static_cast<double>(Lb[b]);
        double K = static_cast<double>(Kb[b]);

        Q += (L / m) - (K / (2.0 * m)) * (K / (2.0 * m));
    }

    return Q;
}

// モジュラリティ (加重)
inline double Modularity_Weighted(const Graph& G,
                                  const std::vector<int>& block_of)
{
    if(block_of.empty())
        return 0.0;

    int nb = 1 + *std::max_element(block_of.begin(), block_of.end());
    if(nb <= 0)
        return 0.0;

    std::vector<double> Kb(nb, 0.0);
    std::vector<double> Lb(nb, 0.0);

    double W = 0.0;

    for_each_undirected_edge(G, [&](int u, int v, double weight) {
        if(u == v)
            return;

        double w = std::abs(weight);
        W += w;
    });

    if(W == 0.0)
        return 0.0;

    for(int u = 0; u < static_cast<int>(num_vertices(G)); ++u){
        int b = block_of[u];

        if(!IsValidBlockIndex(b, nb))
            continue;

        double deg = 0.0;

        for(const auto& e : G.adj[u])
            deg += std::abs(e.weight);

        Kb[b] += deg;
    }

    for_each_undirected_edge(G, [&](int u, int v, double weight) {
        if(u == v)
            return;

        int bu = block_of[u];
        int bv = block_of[v];

        if(IsValidBlockIndex(bu, nb) && bu == bv){
            double w = std::abs(weight);
            Lb[bu] += w;
        }
    });

    double Q = 0.0;

    for(int b = 0; b < nb; ++b){
        double L = Lb[b];
        double K = Kb[b];

        Q += (L / W) - (K / (2.0 * W)) * (K / (2.0 * W));
    }

    return Q;
}

inline std::vector<int>
CountBlocksPerColor(const std::vector<int>& block_color, int nc)
{
    std::vector<int> blocks_per_color(nc, 0);

    for(int c : block_color){
        if(IsValidColorIndex(c, nc))
            blocks_per_color[c] += 1;
    }

    return blocks_per_color;
}

inline std::vector<int>
CountNodesPerColor(const std::vector<int>& nodes_per_block,
                   const std::vector<int>& block_color,
                   int nc)
{
    std::vector<int> nodes_per_color(nc, 0);

    int nb = static_cast<int>(nodes_per_block.size());

    for(int b = 0; b < nb; ++b){
        int c = block_color[b];

        if(IsValidColorIndex(c, nc))
            nodes_per_color[c] += nodes_per_block[b];
    }

    return nodes_per_color;
}

inline BlockEvaluationResult
EvaluatePartitioningMetrics(const Graph& G,
                            const std::vector<int>& block_of,
                            int nb)
{
    BlockEvaluationResult r;

    r.num_vertices = static_cast<int>(num_vertices(G));
    r.num_edges = static_cast<int>(num_edges(G));
    r.num_blocks = nb;

    if(nb <= 0)
        return r;

    r.nodes_per_block = CountNodesPerBlock(block_of, nb);
    r.internal_edges_per_block = CountInternalEdges(G, block_of, nb);

    CountEdgeLocality(G, block_of, nb,
                      r.internal_edges, r.external_edges);

    r.internal_ratio = InternalEdgeRatio(r.internal_edges, r.external_edges);
    r.locality_score = LocalityScore(r.internal_edges, r.external_edges);

    r.avg_nodes_per_block = Mean(r.nodes_per_block);
    r.max_nodes_per_block = *std::max_element(r.nodes_per_block.begin(),
                                              r.nodes_per_block.end());
    r.min_nodes_per_block = *std::min_element(r.nodes_per_block.begin(),
                                              r.nodes_per_block.end());
    r.block_size_stddev = StddevPopulation(r.nodes_per_block);
    r.block_size_gini = GiniCoefficient(r.nodes_per_block);

    r.mean_internal_avg_degree =
        MeanInternalAverageDegree(r.nodes_per_block,
                                  r.internal_edges_per_block);

    r.avg_block_density =
        AverageBlockDensity(r.nodes_per_block,
                            r.internal_edges_per_block);

    Graph T = BuildBlockGraph(G, block_of, BlockEdgeWeight::Binary);

    r.block_degrees = BlockDegreesBinary(T);
    r.block_graph_edges = static_cast<int>(num_edges(T));
    r.block_graph_density = BlockGraphDensity(T);
    r.avg_block_graph_degree = Mean(r.block_degrees);

    r.modularity_unweighted = Modularity_Unweighted(G, block_of);
    r.modularity_weighted = Modularity_Weighted(G, block_of);

    return r;
}

inline BlockEvaluationResult
EvaluatePartitioningMetrics(const Graph& G,
                            const std::vector<int>& block_of,
                            const std::vector<int>& block_color,
                            int nb,
                            int nc)
{
    BlockEvaluationResult r = EvaluatePartitioningMetrics(G, block_of, nb);

    r.num_colors = nc;

    if(nc > 0 && static_cast<int>(block_color.size()) >= nb){
        r.blocks_per_color = CountBlocksPerColor(block_color, nc);
        r.nodes_per_color = CountNodesPerColor(r.nodes_per_block,
                                               block_color,
                                               nc);
    }

    return r;
}

inline void PrintColorStatistics(const BlockEvaluationResult& r)
{
    if(r.num_colors <= 0 || r.blocks_per_color.empty())
        return;

    std::cout << "\n6. Color Statistics" << std::endl;
    std::cout << "   Color"
              << "   Blocks"
              << "   Nodes"
              << "   AvgNodesPerBlock"
              << std::endl;

    for(int c = 0; c < r.num_colors; ++c){
        double avg_nodes = 0.0;

        if(r.blocks_per_color[c] > 0){
            avg_nodes =
                static_cast<double>(r.nodes_per_color[c])
                / r.blocks_per_color[c];
        }

        std::cout << "   "
                  << std::setw(5) << c + 1
                  << "   "
                  << std::setw(6) << r.blocks_per_color[c]
                  << "   "
                  << std::setw(5) << r.nodes_per_color[c]
                  << "   "
                  << std::setw(16)
                  << std::fixed << std::setprecision(2)
                  << avg_nodes
                  << std::endl;
    }
}

inline void PrintPartitioningEvaluation(const BlockEvaluationResult& r)
{
    std::cout << "\n==============================================" << std::endl;
    std::cout << "  Partitioning Evaluation Report" << std::endl;
    std::cout << "==============================================" << std::endl;

    std::cout << "1. Total Blocks          : "
              << r.num_blocks << std::endl;

    std::cout << "   Total Colors          : "
              << r.num_colors << std::endl;

    std::cout << "2. Nodes per Block" << std::endl;

    std::cout << "   - Average             : "
              << std::fixed << std::setprecision(2)
              << r.avg_nodes_per_block << std::endl;

    std::cout << "   - Maximum             : "
              << r.max_nodes_per_block << std::endl;

    std::cout << "   - Minimum             : "
              << r.min_nodes_per_block << std::endl;

    std::cout << "   - Stddev              : "
              << std::setprecision(6)
              << r.block_size_stddev << std::endl;

    std::cout << "   - Gini coefficient    : "
              << std::setprecision(6)
              << r.block_size_gini << std::endl;

    std::cout << "3. Internal Avg Degree   : "
              << std::setprecision(4)
              << r.mean_internal_avg_degree
              << " (Mean of per-block avg)"
              << std::endl;

    std::cout << "4. Block Graph Avg Deg   : "
              << r.avg_block_graph_degree << std::endl;

    std::cout << "5. Modularity" << std::endl;

    std::cout << "   - Unweighted Q        : "
              << std::setprecision(6)
              << r.modularity_unweighted << std::endl;

    std::cout << "   - Weighted Q          : "
              << std::setprecision(6)
              << r.modularity_weighted << std::endl;

    std::cout << "==============================================" << std::endl;

    std::cout << "Internal edges : " << r.internal_edges << "\n";
    std::cout << "External edges : " << r.external_edges << "\n";
    std::cout << "Internal ratio : " << r.internal_ratio << "\n";
    std::cout << "Avg block density : " << r.avg_block_density << "\n";
    std::cout << "Block graph edges : " << r.block_graph_edges << "\n";
    std::cout << "Block graph density : " << r.block_graph_density << "\n";
    std::cout << "Locality score : " << r.locality_score << "\n";
    // std::cout << "Block size stddev : " << r.block_size_stddev << "\n";
    // std::cout << "Block size Gini   : " << r.block_size_gini << "\n";

    PrintColorStatistics(r);
}

/**
 * 旧 eblock.cpp と同じ呼び出し口。
 *
 * 今後の Leiden/CPM 反復プログラムでは、表示が不要な場合は
 * EvaluatePartitioningMetrics() を直接呼び、BlockEvaluationResult の値を
 * 停止条件や gamma 更新条件に使うとよい。
 */
inline void EvaluatePartitioning(const Graph& G,
                                 const std::vector<int>& block_of,
                                 const std::vector<int>& block_color,
                                 int nb,
                                 int nc)
{
    BlockEvaluationResult r =
        EvaluatePartitioningMetrics(G, block_of, block_color, nb, nc);

    PrintPartitioningEvaluation(r);
}

// 色情報がまだ無い段階の評価用。
inline void EvaluatePartitioning(const Graph& G,
                                 const std::vector<int>& block_of,
                                 int nb)
{
    BlockEvaluationResult r =
        EvaluatePartitioningMetrics(G, block_of, nb);

    PrintPartitioningEvaluation(r);
}

// 旧関数名を残した薄いラッパー群。
inline void EvaluateEdgeLocality(const Graph& G,
                                 const std::vector<int>& block_of,
                                 int nb)
{
    int internal = 0;
    int external = 0;

    CountEdgeLocality(G, block_of, nb, internal, external);

    std::cout << "Internal edges : " << internal << "\n";
    std::cout << "External edges : " << external << "\n";
    std::cout << "Internal ratio : "
              << InternalEdgeRatio(internal, external) << "\n";
}

inline void EvaluateBlockDensity(const Graph& G,
                                 const std::vector<int>& block_of,
                                 int nb)
{
    std::vector<int> nodes_per_block = CountNodesPerBlock(block_of, nb);
    std::vector<int> internal_edges = CountInternalEdges(G, block_of, nb);

    std::cout << "Avg block density : "
              << AverageBlockDensity(nodes_per_block, internal_edges)
              << "\n";
}

inline void EvaluateBlockGraph(const Graph& T)
{
    std::cout << "Block graph edges : "
              << num_edges(T) << "\n";

    std::cout << "Block graph density : "
              << BlockGraphDensity(T) << "\n";
}

// inline void EvaluateLocalityScore(const Graph& G,
//                                   const std::vector<int>& block_of)
// {
//     int nb = 1 + *std::max_element(block_of.begin(), block_of.end());
//     int internal = 0;
//     int external = 0;

//     CountEdgeLocality(G, block_of, nb, internal, external);

//     std::cout << "Locality score : "
//               << LocalityScore(internal, external) << "\n";
// }
void EvaluateLocalityScore(
    const Graph& G,
    const std::vector<int>& block_of)
{
    int internal = 0;
    int external = 0;

    int N = static_cast<int>(num_vertices(G));

    for_each_undirected_edge(G, [&](int u, int v, double) {
        if(block_of[u] == block_of[v])
            internal++;
        else
            external++;
    });

    double score =
        static_cast<double>(internal)
        / static_cast<double>(external + 1);

    std::cout
        << "Locality score : "
        << score << "\n";
}

inline void EvaluateLoadBalance(const std::vector<int>& nodes_per_block)
{
    std::cout << "Block size stddev : "
              << StddevPopulation(nodes_per_block) << "\n";

    std::cout << "Block size Gini   : "
              << GiniCoefficient(nodes_per_block) << "\n";
}

inline void EvaluateColorStatistics(const std::vector<int>& nodes_per_block,
                                    const std::vector<int>& block_color,
                                    int nc)
{
    BlockEvaluationResult r;

    r.num_blocks = static_cast<int>(nodes_per_block.size());
    r.num_colors = nc;
    r.nodes_per_block = nodes_per_block;
    r.blocks_per_color = CountBlocksPerColor(block_color, nc);
    r.nodes_per_color = CountNodesPerColor(nodes_per_block, block_color, nc);

    PrintColorStatistics(r);
}
