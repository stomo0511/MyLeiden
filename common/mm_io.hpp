#pragma once
#include <string>
#include <vector>

#include "Types.hpp"

Graph Read_MM_UD(const std::string& file_name);
void WriteReorderedMatrixMarket_UD(
    const Graph& G,
    const std::vector<int>& perm,
    const std::string& file_name);
void WritePermutation_1Based(
    const std::vector<int>& perm,
    const std::string& file_name);
std::string file_stem(const std::string& path);