#include "Leiden.hpp"
#include "common/BlockIO.hpp"
#include "common/BlockTypes.hpp"
#include "common/Coloring.hpp"
#include "common/MM_IO.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

#if defined(LEIDEN_DRIVER_CPM)
constexpr const char* kMethodName = "LeidenCP";
constexpr const char* kResolutionExample = "0.001";
#elif defined(LEIDEN_DRIVER_MODULARITY)
constexpr const char* kMethodName = "LeidenMD";
constexpr const char* kResolutionExample = "1.0";
#else
#error "Define either LEIDEN_DRIVER_CPM or LEIDEN_DRIVER_MODULARITY"
#endif

void PrintUsage(const char* program_name)
{
    std::cerr << "Usage:\n"
              << "  " << program_name << " <matrix.mm> <resolution>\n\n"
              << "Example:\n"
              << "  " << program_name << " Saad.mm "
              << kResolutionExample << "\n";
}

bool FileExists(const std::string& path)
{
    std::ifstream in(path);
    return static_cast<bool>(in);
}

bool IsSafeFilenameToken(const std::string& token)
{
    if (token.empty()) {
        return false;
    }
    for (char ch : token) {
        const bool ok =
            (ch >= '0' && ch <= '9') ||
            (ch >= 'a' && ch <= 'z') ||
            (ch >= 'A' && ch <= 'Z') ||
            ch == '.' || ch == '_' || ch == '-' || ch == '+';
        if (!ok) {
            return false;
        }
    }
    return true;
}

double ParsePositiveResolution(const std::string& text)
{
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("resolution is not a valid number: " + text);
    }
    if (consumed != text.size()) {
        throw std::invalid_argument("resolution has trailing characters: " + text);
    }
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument("resolution must be finite and positive");
    }
    return value;
}

void SetBinaryEdgeWeights(Graph& G)
{
    for (std::vector<Edge>& adjacency : G.adj) {
        for (Edge& e : adjacency) {
            e.weight = 1.0;
        }
    }
}

std::unique_ptr<QualityFunction> MakeQualityFunction(double resolution)
{
#if defined(LEIDEN_DRIVER_CPM)
    return std::unique_ptr<QualityFunction>(
        new CPMQualityFunction(resolution));
#elif defined(LEIDEN_DRIVER_MODULARITY)
    return std::unique_ptr<QualityFunction>(
        new ModularityQualityFunction(resolution));
#endif
}

BlockPartition MakeBlockPartitionFromLeiden(const LeidenPartition& partition,
                                            int n)
{
    if (static_cast<int>(partition.community_of.size()) != n) {
        throw std::invalid_argument("partition size does not match graph");
    }

    BlockPartition block_partition;
    block_partition.n = n;
    block_partition.s = 1;
    block_partition.block_of.assign(n, -1);

    int max_block = -1;
    for (int v = 0; v < n; ++v) {
        const int block = partition.community_of[v];
        if (block < 0) {
            throw std::runtime_error("Leiden returned a negative community id");
        }
        block_partition.block_of[v] = block;
        max_block = std::max(max_block, block);
    }

    block_partition.nb = max_block + 1;
    block_partition.blocks.assign(block_partition.nb, {});
    for (int v = 0; v < n; ++v) {
        const int block = block_partition.block_of[v];
        if (block >= block_partition.nb) {
            throw std::runtime_error("block id is out of range");
        }
        block_partition.blocks[block].push_back(v);
    }

    for (int block = 0; block < block_partition.nb; ++block) {
        if (block_partition.blocks[block].empty()) {
            throw std::runtime_error("Leiden final partition is not compact");
        }
    }
    return block_partition;
}

void ValidateColoring(const Graph& block_graph,
                      const std::vector<int>& block_color,
                      int num_colors)
{
    if (num_colors <= 0 && num_vertices(block_graph) > 0) {
        throw std::runtime_error("coloring produced no colors");
    }
    if (static_cast<int>(block_color.size()) != num_vertices(block_graph)) {
        throw std::runtime_error("color vector size does not match block graph");
    }

    for (int block = 0; block < num_vertices(block_graph); ++block) {
        const int color = block_color[block];
        if (color < 0 || color >= num_colors) {
            throw std::runtime_error("block color is out of range");
        }
    }

    for_each_undirected_edge(block_graph, [&](int u, int v, double) {
        if (u != v && block_color[u] == block_color[v]) {
            throw std::runtime_error("invalid coloring: adjacent blocks share a color");
        }
    });
}

} // namespace

int main(int argc, char** argv)
{
    const char* program_name = (argc > 0) ? argv[0] : kMethodName;
    if (argc != 3) {
        PrintUsage(program_name);
        return EXIT_FAILURE;
    }

    const std::string matrix_file = argv[1];
    const std::string resolution_text = argv[2];

    try {
        if (!FileExists(matrix_file)) {
            throw std::runtime_error("input file does not exist: " + matrix_file);
        }
        if (!IsSafeFilenameToken(resolution_text)) {
            throw std::invalid_argument(
                "resolution contains characters that are unsafe for output filenames");
        }

        const double resolution = ParsePositiveResolution(resolution_text);
        Graph G = Read_MM_UD(matrix_file);
        SetBinaryEdgeWeights(G);
        const LeidenGraphStats stats = BuildLeidenGraphStats(G);
        const std::unique_ptr<QualityFunction> quality =
            MakeQualityFunction(resolution);

        LeidenOptions options;
        options.theta = 0.01;
        options.seed = 0;
        options.max_levels = 0;

        const LeidenResult result = Leiden(G, stats, *quality, options);
        const BlockPartition block_partition =
            MakeBlockPartitionFromLeiden(result.partition, num_vertices(G));

        const Graph block_graph =
            BuildBlockGraph(G,
                            block_partition.block_of,
                            BlockEdgeWeight::Binary);
        std::vector<int> block_color;
        int num_colors = Greedy_Coloring(block_graph, block_color);
        RelabelColorsByClassSize(block_color);
        num_colors =
            block_color.empty()
                ? 0
                : 1 + *std::max_element(block_color.begin(),
                                        block_color.end());
        ValidateColoring(block_graph, block_color, num_colors);

        const std::string base = file_stem(matrix_file) + "_" +
                                 kMethodName + "_gamma" + resolution_text;
        const std::string block_file = base + ".blk";
        const std::string color_file = base + ".bcol";

        WriteBlockInfo_1Based(block_partition.block_of, block_file);
        WriteBlockColor_1Based(block_color, num_colors, color_file);

        const double final_quality =
            quality->quality(G, stats, result.partition);

        std::cout << "Method              : " << kMethodName << "\n"
                  << "Matrix              : " << matrix_file << "\n"
                  << "Resolution          : " << resolution_text << "\n"
                  << "Theta               : " << options.theta << "\n"
                  << "Seed                : " << options.seed << "\n"
                  << "Number of vertices  : " << num_vertices(G) << "\n"
                  << "Number of edges     : " << num_edges(G) << "\n"
                  << "Number of blocks    : " << block_partition.nb << "\n"
                  << "Number of colors    : " << num_colors << "\n"
                  << "Leiden levels       : " << result.num_levels << "\n"
                  << "Leiden moves        : " << result.total_moves << "\n"
                  << "Quality value (internal objective): "
                  << final_quality << "\n"
                  << "Block file          : " << block_file << "\n"
                  << "Color file          : " << color_file << "\n";

        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        PrintUsage(program_name);
        return EXIT_FAILURE;
    }
}
