#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstddef>
#include "Types.hpp"

// 頂点と次数のペアを格納する構造体
struct VertexDegree {
    Vertex v;
    std::size_t degree;

    bool operator<(const VertexDegree& other) const {
        return degree > other.degree;  // 降順（largest first）
    }
};

enum class BlockPolicy {
    FIFO = 1,        // Policy 1: C をFIFOで取り出す
    MaxEdges = 2,    // Policy 2: ブロック内への辺本数が最大
    Weighted = 3     // Policy 3: ブロック内ノードとの|a_ij|合計が最大
};

static const char* policy_name(BlockPolicy policy) {
    switch (policy) {
        case BlockPolicy::FIFO:      return "Policy 1 (FIFO)";
        case BlockPolicy::MaxEdges:  return "Policy 2 (MaxEdges)";
        case BlockPolicy::Weighted:  return "Policy 3 (Weighted)";
        default:                     return "Unknown Policy";
    }
};

struct BlockPartition {
    int n = 0;                           // 頂点数
    int s = 1;                           // 目標ブロックサイズ
    int nb = 0;                          // 生成ブロック数
    std::vector<int> block_of;           // 各ノード→ブロックID (-1: 未割当)
    std::vector<std::vector<int>> blocks; // 各ブロックのノード列（元番号のまま）
};

Graph Read_MM_UD(const std::string& file_name);
void WriteReorderedMatrixMarket_UD(
    const Graph& G,
    const std::vector<int>& perm,
    const std::string& file_name);
void WritePermutation_1Based(
    const std::vector<int>& perm,
    const std::string& file_name);
std::string file_stem(const std::string& path);