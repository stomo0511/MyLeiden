#include "mm_io.hpp"

#include <algorithm>

Graph Read_MM_UD(const std::string& file_name)
{
    std::ifstream infile(file_name);
    if (!infile)
    {
        std::cerr << "Error: Cannot open file " << file_name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    // コメント行をスキップし、サイズを読み取る
    do
    {
        getline(infile, line);
    } while (!infile.eof() && line[0] == '%');

    std::istringstream size_stream(line);
    int row_size, col_size, num_elements;
    size_stream >> row_size >> col_size >> num_elements;

    if (row_size != col_size)
    {
        std::cerr << "Error: The matrix is not square." << std::endl;
        exit(EXIT_FAILURE);
    }

    Graph G(row_size);

    int row, col;
    double val;
    while (getline(infile, line))
    {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream edge_stream(line);
        val = 1.0;
        edge_stream >> row >> col;
        if (!(edge_stream >> val)) val = 1.0;

        // Matrix Market is 1-based index → convert to 0-based
        if (row > col)
        {
            add_undirected_edge(G, row - 1, col - 1, val);
        }
    }

    return G;
}

void WriteReorderedMatrixMarket_UD(
    const Graph& G,
    const std::vector<int>& perm,
    const std::string& file_name)
{
    const int n = G.n;

    // 無向グラフとして、upper triangle のみ出力
    struct Entry {
        int i, j;
        double a;
    };
    std::vector<Entry> entries;
    entries.reserve(num_edges(G));

    for_each_undirected_edge(G, [&](int u, int v, double w) {
        int ru = perm[u];
        int rv = perm[v];
        if (ru > rv) std::swap(ru, rv);

        // Matrix Market は 1-based
        entries.push_back({ru + 1, rv + 1, w});
    });

    std::sort(entries.begin(), entries.end(),
              [](const Entry& x, const Entry& y) {
                  if (x.i != y.i) return x.i < y.i;
                  return x.j < y.j;
              });

    std::ofstream ofs(file_name);
    if (!ofs) {
        throw std::runtime_error("failed to open output file: " + file_name);
    }

    ofs << "%%MatrixMarket matrix coordinate real symmetric\n";
    ofs << n << " " << n << " " << entries.size() << "\n";
    for (const auto& e : entries) {
        ofs << e.i << " " << e.j << " " << e.a << "\n";
    }
}

// 順列をテキストでも保存しておくと後で便利
// 1行: old_index(1-based) new_index(1-based)
void WritePermutation_1Based(
    const std::vector<int>& perm,
    const std::string& file_name)
{
    std::ofstream ofs(file_name);
    if (!ofs) {
        throw std::runtime_error("failed to open permutation file: " + file_name);
    }

    for (int old_id = 0; old_id < static_cast<int>(perm.size()); ++old_id) {
        ofs << (old_id + 1) << " " << (perm[old_id] + 1) << "\n";
    }
}

// ---- ユーティリティ：入力パスから拡張子を外したstemを得る
std::string file_stem(const std::string& path) {
    // 末尾の '/' or '\' の次からファイル名部分
    size_t pos = path.find_last_of("/\\");
    std::string name = (pos == std::string::npos) ? path : path.substr(pos+1);
    // 拡張子を落とす
    size_t dot = name.find_last_of('.');
    return (dot == std::string::npos) ? name : name.substr(0, dot);
}
