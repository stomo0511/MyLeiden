#pragma once

#include <cstddef>
#include <vector>

struct Edge {
    int to;
    double weight;
};

struct SimpleGraph {
    int n = 0;
    std::vector<std::vector<Edge>> adj;

    SimpleGraph() = default;
    explicit SimpleGraph(int n_) : n(n_), adj(n_) {}
};

using Graph = SimpleGraph;
using Vertex = int;

// Block graph T を作る（binary/count/abs-sum を選べる）
enum class BlockEdgeWeight { Binary, Count, AbsSum };

inline Graph MakeGraph(int n)
{
    return Graph(n);
}

inline int num_vertices(const Graph& G)
{
    return G.n;
}

inline std::size_t degree(int u, const Graph& G)
{
    return G.adj[u].size();
}

inline void add_undirected_edge(Graph& G, int u, int v, double weight = 1.0)
{
    G.adj[u].push_back({v, weight});
    if (u != v)
        G.adj[v].push_back({u, weight});
}

inline std::size_t num_edges(const Graph& G)
{
    std::size_t m = 0;
    for (const auto& nbrs : G.adj)
        m += nbrs.size();
    return m / 2;
}

inline bool find_edge_weight(const Graph& G, int u, int v, double& weight)
{
    for (const auto& e : G.adj[u]) {
        if (e.to == v) {
            weight = e.weight;
            return true;
        }
    }
    return false;
}

template <class Func>
inline void for_each_undirected_edge(const Graph& G, Func&& func)
{
    for (int u = 0; u < G.n; ++u) {
        for (const auto& e : G.adj[u]) {
            if (u <= e.to)
                func(u, e.to, e.weight);
        }
    }
}
