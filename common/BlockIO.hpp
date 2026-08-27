#pragma once
#include <iostream>
#include <fstream>

void WriteBlockInfo_1Based(const std::vector<int>& block_of, const std::string& out_path);
void WriteBlockColor_1Based(const std::vector<int>& block_color, int nc, const std::string& out_path);
Graph BuildBlockGraph(const Graph& G,
                      const std::vector<int>& block_of,
                      BlockEdgeWeight mode);

int ReadBlockInfo_1Based(const std::string& path, int N, std::vector<int>& block_of);
int ReadBlockColor_1Based(const std::string& path, int nb, std::vector<int>& block_color); 
