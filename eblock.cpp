#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "common/Types.hpp"
#include "common/MM_IO.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

int main(int argc, char** argv)
{
    if(argc < 4){
        std::cerr
            << "Usage: "
            << argv[0]
            << " <matrix.mtx> <block_file.blk> <color_file.bcol>"
            << std::endl;

        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);
    int N = num_vertices(G);

    std::vector<int> block_of;
    int nb = ReadBlockInfo_1Based(argv[2], N, block_of);

    std::vector<int> block_color;
    int nc = ReadBlockColor_1Based(argv[3], nb, block_color);

    EvaluatePartitioning(G, block_of, block_color, nb, nc);

    /*
     * Additional load imbalance metrics:
     *
     * 1. Maximum block size / average block size
     * 2. Maximum color-wise node count / average color-wise node count
     *
     * Assumption:
     *   ReadBlockInfo_1Based() and ReadBlockColor_1Based() read 1-based
     *   files and convert the indices to 0-based internal representation.
     */

    std::vector<int> block_size(nb, 0);

    for(int v = 0; v < N; ++v){
        int b = block_of[v];

        if(b < 0 || b >= nb){
            std::cerr
                << "Error: invalid block index "
                << b
                << " for vertex "
                << v
                << std::endl;
            return 1;
        }

        block_size[b]++;
    }

    int max_block_size = 0;
    if(nb > 0){
        max_block_size = *std::max_element(block_size.begin(), block_size.end());
    }

    double avg_block_size = 0.0;
    if(nb > 0){
        avg_block_size = static_cast<double>(N) / static_cast<double>(nb);
    }

    double max_block_over_avg = 0.0;
    if(avg_block_size > 0.0){
        max_block_over_avg =
            static_cast<double>(max_block_size) / avg_block_size;
    }

    std::vector<int> color_node_count(nc, 0);

    for(int b = 0; b < nb; ++b){
        int c = block_color[b];

        if(c < 0 || c >= nc){
            std::cerr
                << "Error: invalid color index "
                << c
                << " for block "
                << b
                << std::endl;
            return 1;
        }

        color_node_count[c] += block_size[b];
    }

    int max_color_node_count = 0;
    if(nc > 0){
        max_color_node_count =
            *std::max_element(color_node_count.begin(), color_node_count.end());
    }

    double avg_color_node_count = 0.0;
    if(nc > 0){
        avg_color_node_count = static_cast<double>(N) / static_cast<double>(nc);
    }

    double max_color_over_avg = 0.0;
    if(avg_color_node_count > 0.0){
        max_color_over_avg =
            static_cast<double>(max_color_node_count) / avg_color_node_count;
    }

    std::cout << std::fixed << std::setprecision(6);

    std::cout << std::endl;
    std::cout << "Additional Load Imbalance Metrics" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Max block size / Avg block size   : "
              << max_block_over_avg << std::endl;
    std::cout << "Max color nodes / Avg color nodes : "
              << max_color_over_avg << std::endl;

    return 0;
}