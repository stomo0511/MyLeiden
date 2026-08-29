# MyLeiden 並列化・最適化 開発履歴 / Cursor 引き継ぎノート

## 0. この文書の目的

この文書は、`MyLeiden` プロジェクトで進めてきた Leiden 法の C++ 実装、逐次最適化、OpenMP 並列化、Stage-4A / Stage-4B 系の検証・性能評価の経緯を、Cursor 上で継続開発するためにまとめた引き継ぎ資料である。

Cursor ではモデルとして **GPT-5.6 Sol** を使用する想定である。

今後の会話・実装では、この文書をプロジェクトの背景情報として参照し、過去に検証済みの事項を不用意にやり直さないこと。

---

# 1. プロジェクトの目的

`MyLeiden` は、疎行列の非零パターンから作成した隣接グラフに Leiden 法を適用し、得られたコミュニティをブロックとして利用する C++ 実装である。

主な用途は、疎行列ソルバ、特に ICCG 法の並列前処理におけるブロック生成である。

主な比較対象は以下。

- ABMC
- Leiden + Modularity (`LeidenMD`)
- Leiden + CPM (`LeidenCP`)

Leiden で生成した block partition に対して greedy coloring を行い、block multi-coloring 型の並列前処理へ利用する。

研究上の重要な観点は、Leiden 単体の速度だけではなく、partition quality、block 数、block color 数、block size imbalance、ICCG iteration count、cache locality、solver execution time まで含めた downstream performance である。

---

# 2. 基準アルゴリズム

参照している Leiden 法は Traag, Waltman, van Eck, “From Louvain to Leiden: guaranteeing well-connected communities”, Scientific Reports, 2019。Supplement Algorithm A.2 を基準としている。

概略:

```text
Leiden
  P <- MoveNodesFast
  if P is not done:
      Prefined <- RefinePartition
      G <- AggregateGraph(G, Prefined)
      P <- coarse partition based on original non-refined P
  repeat
  flatten to original vertices
```

主要フェーズ:

1. `MoveNodesFast`
2. `RefinePartition`
3. `AggregateGraph`
4. coarse partition construction
5. multilevel loop
6. flatten

---

# 3. 主なソース構成

```text
MyLeiden/
  Makefile
  Leiden.hpp
  Leiden.cpp
  LeidenTypes.hpp
  QualityFunction.hpp
  QualityFunction.cpp
  leiden_main.cpp
  eblock.cpp

  common/
    Types.hpp
    BlockTypes.hpp
    MM_IO.hpp
    MM_IO.cpp
    BlockIO.hpp
    BlockIO.cpp
    Coloring.hpp
    Coloring.cpp
    Block_Eval.hpp
    Timer.hpp

  tests/
    test_quality.cpp
    test_move_nodes_fast.cpp
    test_refinement.cpp
    test_aggregate.cpp
    test_coarse_partition.cpp
    test_leiden.cpp
    test_block_eval.cpp
    test_move_nodes_fast_stage4b_tsan.cpp
    test_openmp_tsan_runtime.cpp
```

現在の開発 branch:

```text
Parallel4
```

---

# 4. Graph 表現

`common/Types.hpp` では独自の adjacency-list graph を使用している。

```cpp
struct Edge {
    int to;
    double weight;
};

struct SimpleGraph {
    int n = 0;
    std::vector<std::vector<Edge>> adj;
};

using Graph = SimpleGraph;
using Vertex = int;
```

重要事項:

- self-loop を扱う。
- parallel edge を許容する。
- `num_edges()` の self-loop 数え上げ問題は修正済み。
- Leiden aggregation では edge pair を集約してから coarse graph を作る。
- `BuildBlockGraph` が `bk == bl` を除外するのは coloring 用として正しい。
- self-loop / parallel-edge corner cases はテスト済み。

---

# 5. QualityFunction abstraction

Leiden の quality function は抽象化されている。主に Modularity と CPM を実装済み。

## CPM move delta

```text
delta =
    (weight_to_target - weight_to_source)
    - gamma * node_size[v]
      * (target_size - source_size_without_v)
```

refinement で source が singleton の場合:

```text
delta =
    weight_to_target
    - gamma * node_size[v] * target_size
```

## Modularity move delta

```text
delta =
    2 * (weight_to_target - weight_to_source)
    - gamma * kv
      * (target_strength - source_strength_without_v)
      / total_edge_weight
```

Stage-4B では asynchronous evaluation 用に scalar snapshot API を追加。

```text
supportsConcurrentMoveEvaluation()
deltaMoveFromCommunitySnapshot(...)
```

CPM / Modularity は既存 `deltaMoveFromWeights()` と数学的に同一の式を使用する。未知の `QualityFunction` は serial path に fallback する。

---

# 6. Stage 1〜5: Leiden 基本実装

以下は実装・テスト済み。

1. Quality functions / graph stats
2. `MoveNodesFast`
3. `RefinePartition`
4. `AggregateGraph`
5. `BuildCoarsePartition`
6. full multilevel Leiden
7. flatten to original graph

テスト対象:

- single vertex
- no-edge graph
- self-loop
- weighted edge
- parallel edge
- CPM
- Modularity
- seed reproducibility
- different seeds
- >= 3 multilevel levels
- max_levels
- empty graph
- quality preservation
- aggregation correctness
- coarse partition correctness

ASan / UBSan も複数段階で確認済み。

---

# 7. Driver

主な実行形:

```bash
./LeidenCP <matrix.mm> <resolution> [--debug] [--debug-interval N]
./LeidenMD <matrix.mm> <resolution> [--debug] [--debug-interval N]
```

処理:

1. Matrix Market 読み込み
2. binary edge 化
3. graph stats
4. Leiden
5. block graph
6. greedy coloring
7. `RelabelColorsByClassSize`
8. coloring validation
9. `.blk` / `.bcol` 出力

Execution time は matrix input と output file write を除き、Leiden + coloring 等を含む。

---

# 8. 逐次版で行った主な最適化

## 8.1 `MakePartition` の O(nc*m) 問題

旧実装では community ごとに全 edge を走査して内部 edge weight を求めていた。これを一度の `for_each_undirected_edge(...)` で集計する形に変更し、実質 `O(n + nc + m)` に改善。

## 8.2 `RefinementCommunityStats`

subset 内の community 情報を `member_count`, `mass`, `external_weight` として保持。subset-local `unordered_map` を利用。

## 8.3 subset membership marker

marker/generation buffer を `RefinePartition` 内で再利用し、allocation を削減。

## 8.4 `NeighborCommunityScratch`

neighbor community weights 構築時の一時 `unordered_map` allocation を避けるため、`weights`, `marks`, `generation`, `touched` を持つ reusable scratch を実装。

## 8.5 `MoveNodeToCommunityFromWeights`

move の際に source / target の neighbor weight を再走査せず、計算済みの `weight_to_source`, `weight_to_target`, `self_loop_weight` を渡して stats 更新する fast path を追加。

## 8.6 Empty community 管理

当初 `std::set<Community>` を使っていたが、profiling で move cost の大きな割合を占めたため array + smallest-empty-ID hint 型へ変更。逐次実装は igraph 相当の速度に近づき、細かな逐次最適化はいったん停止。

---

# 9. AggregateGraph OpenMP 並列化

## Stage 1

global edge list + sort/reduce 型を試したが、1 thread 時の overhead が大きかった。

Emilia_923:

```text
1T ~1.67 s
2T ~1.63 s
4T ~1.57 s
8T ~1.78 s
```

## Stage 2

community-local accumulation / scratch 方式へ変更。

```text
1T ~0.204 s
2T ~0.162 s
4T ~0.089 s
8T ~0.071 s
```

AggregateGraph 並列化は完了扱い。

---

# 10. RefinePartition OpenMP 並列化: Stage 3

parent community ごとの refinement subset を OpenMP で並列化。static と `dynamic,1` を比較。

Emilia:

```text
threads   static      dynamic,1
1         9.11 s      9.33 s
2         6.05 s      5.44 s
4         4.46 s      4.76 s
8         5.23 s      4.32 s
```

production default は static。Stage-3 では thread count 間で partition output が byte-identical になる deterministic subset RNG を導入。

---

# 11. Stage-4A: Parallel proposal + serial deterministic commit

`MoveNodesFast` 並列化の最初の試み。

```text
affected_current
parallel immutable-snapshot proposal generation
serial deterministic commit
revalidation
neighbor reactivation
affected_next
repeat rounds
```

特徴:

- parallel proposal 中は partition immutable
- commit は逐次
- stale proposal を revalidate
- affected set を double-buffer
- deterministic order で commit

race-free を優先した保守的設計だったが、round-based work amplification、repeated proposals、serial commit、revalidation、full-n scan が重かった。

---

# 12. Stage-4A.1: affected-set / reactivation 改善

serial `MoveNodesFast` の queue semantics に近づけるため reactivation rule を改善。

move `v -> B` 後:

- `v` 自身は自己移動だけでは再activateしない
- non-self neighbor `u` について `community_of[u] != B` のときだけ activate
- rejected/stale proposal は自己 reactivation
- no-positive proposal は自己 reactivation しない
- empty-target collision は reactivation

Emilia_923 8T 代表値:

```text
visits              12,377,795
moves                 2,563,864
rounds                       601
positive proposals     2,897,436
rejected                 333,572
MoveNodesFast           ~13.63 s
```

Stage-3 serial:

```text
visits 3,742,978
moves  1,307,167
MoveNodesFast ~5.7 s （当時）
```

work amplification:

```text
visits ~3.31x
moves  ~1.96x
```

frontier compaction も試したが悪化したため reject。

---

# 13. RefinePartition 巨大 slowdown の原因解析

Stage-4A 系 partition では `RefinePartition` が異常に遅くなるケースが発生。

主要因:

- parent subset 内 active community 数増大
- active-community full scan
- candidate processing
- 巨大 subset critical path

重要でなかったもの:

- subset count
- refinement move count
- scratch init
- MakePartition
- 単純な adjacency volume 増加だけ

active-community iteration と delta evaluation が約100倍規模まで増えたケースがあった。

---

# 14. Exact sparse refinement targets

refinement candidate target の full active-community scan を削減。

現在の CPM / Modularity の式では、source singleton refinement 時に正の mass を持つ non-neighbor target は delta < 0 になるため、exact candidate set を次に限定可能。

```text
neighbor communities
+ source
+ active nonpositive-mass exceptional communities
```

zero mass / zero total edge weight / unknown QualityFunction などは fallback/full scan を維持。candidate order は sort + unique で deterministic にし、`discrete_distribution` の RNG semantics を維持。

Emilia_923:

```text
full active iterations      ~4.27 billion
sparse target iterations    ~17.4 million
full delta evaluations      ~4.27 billion
sparse delta evaluations    ~16.5 million
```

iteration 数は約250倍減ったが wall-time 改善は約1.6倍程度で、別 bottleneck が残っていた。

---

# 15. Incremental active-list optimization

profiling により、refinement move ごとの `RefreshActiveCommunities(...)` が最大級の bottleneck と判明。

Emilia level 1:

```text
move refresh calls       872,025
entries scanned        8.16 billion
active entries output  4.27 billion
refresh CPU sum        ~151 s
```

`RefinementCommunityEntry` に位置情報を追加。

```cpp
std::size_t active_position;
std::size_t nonpositive_mass_position;
```

active list を swap-remove + push_back で incrementally 更新。

```text
source deactivates iff old_source_count == 1
target activates iff old_target_count == 0
```

zero/nonpositive-mass list も source/target の状態変化だけ更新。

比較用 macro:

```text
REFINE_FULL_REFRESH_ACTIVE_MAINTENANCE
REFINE_ZERO_MASS_ACTIVE_ONLY
```

Emilia_923:

```text
RefinePartition:
42.98 s -> 1.774 s

largest subset:
13.27 s -> 0.109 s

incremental update CPU:
~0.072 s
```

約24倍改善。exact sparse + incremental active-list の組み合わせで旧 Refine bottleneck は解消。

---

# 16. Stage-4A.1 + incremental Refine の代表結果

後日の同一 GCC/O3 環境で再測定した Emilia_923 Stage-4A.1 8T median:

```text
MoveNodesFast     12.724 s
RefinePartition    1.108 s
Leiden total      14.617 s
processed         12,377,795
moves              2,563,864
quality           3.60897e7
blocks/colors     33/8
```

現時点の production default は Stage-4A.1。

---

# 17. Stage-4B: asynchronous race-free parallel MoveNodesFast

Stage-4A の round + serial commit architecture を捨て、非同期 local moving を実験実装。

設計:

```text
atomic community labels
community occupancy
per-community locks
atomic affected/queued state
shared work queue
outstanding counter
thread-local neighbor scratch
optimistic candidate evaluation
source/target locked revalidation
immediate commit
neighbor reactivation
verification sweeps
final partition rebuild
```

global lock は使用しない。

---

# 18. ConcurrentMoveState

| State | 管理 |
|---|---|
| `community_of[v]` | atomic |
| occupancy | atomic、commit lock 内で再検証 |
| current community mass | community mutex |
| source/target locks | community ID 昇順 |
| affected/queued state | atomic 3-state |
| work queue | short mutex |
| outstanding work | atomic |
| neighbor weights | thread-local |
| candidate list | thread-local |
| internal edge weight | 終了後 rebuild |
| community size/strength | 終了後 rebuild |
| empty flags/hint | occupancyから管理、終了後 rebuild |

source/target lock は `min(source,target)` → `max(source,target)` の順で取得し deadlock を防止。

---

# 19. Stage-4B queue state machine

```text
0 = IDLE
1 = QUEUED or PROCESSING
2 = REQUEUE_REQUESTED
```

```text
outstanding = physically queued vertices + processing vertices
```

enqueue:

```text
IDLE -> QUEUED CAS
outstanding++
queue lock
queue.push_back()
```

processing completion:

```text
previous = state.exchange(IDLE)
if previous == REQUEUE_REQUESTED:
    enqueue(v)
outstanding--
```

phase 終了時に以下を runtime check。

```text
queue empty
outstanding == 0
all vertex states == IDLE
```

---

# 20. Stage-4B.1: Race Isolation and TSan Validation

当初 GCC TSan + libgomp で多数の warning が出たため性能測定前に停止。

Stage-4B 専用 TSan executable を作り、Stage-4A / Refine / Aggregate を除外。さらに OpenMP と `std::thread` の両 backend から同一 worker logic を実行する形へ整理。

結果:

```text
std::thread + TSan:
1T clean
2T clean
4T clean
8T clean

OpenMP/libgomp:
warning remains
```

warning は主に queue/shared stack object lifetime、outstanding atomic lifetime、worker-local vector lifetime、region 終了後の profile/rebuild access、shared object destruction に集中。

core candidate/commit/requeue logic の unsynchronized plain access は確認されなかった。

---

# 21. Stage-4B.2: GCC/libgomp TSan interoperability validation

環境:

```text
Ubuntu 22.04.5 LTS
Kernel 5.15.0-186-generic
GCC/G++ 11.4.0
libgomp.so.1
libgomp package 12.3.0-1ubuntu1~22.04.3
libtsan.so.0
```

TSan matrix:

```text
std::thread 1T: 0 warning
std::thread 2T: 0
std::thread 4T: 0
std::thread 8T: 0

OpenMP 1T: 0
OpenMP 2T: ~18
OpenMP 4T: ~19
OpenMP 8T: ~16
```

## Stage-4B 非依存 minimal OpenMP program

安全な最小プログラムで同型の warning を再現。

```text
1T 0 warnings
2T 6 warnings
4T 6 warnings
8T 6 warnings
```

よって GCC/libgomp–TSan interoperability issue と強く整合する。ただし「confirmed false positive」とは断定しない。

Stage-4B.2 結論:

```text
Stage-4B core worker race-free evidence: strong
OpenMP backend race-free evidence: strong
Stage-4B correctness: verified
Performance benchmark: ready
```

---

# 22. Stage-4B.3: Performance Evaluation

環境:

```text
Host: epyc
OS: Ubuntu 22.04.5
CPU: AMD EPYC 7543 32-Core Processor
1 socket
32 cores
1 thread/core
NUMA nodes: 1
L3: 256 MiB

GCC 11.4.0
-O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp
OMP_DYNAMIC=FALSE
OMP_PROC_BIND unset
OMP_PLACES unset
```

warm-up 1回を除外し各条件5回、primary timing は median。

---

# 23. Emilia_923: Stage-4B timing

`LeidenMD`, resolution = 1.0。

| Threads | Move median | Leiden median |
|---:|---:|---:|
| 1 | 10.569 s | 19.575 s |
| 2 | 9.325 s | 44.108 s |
| 4 | 6.999 s | 42.240 s |
| 8 | 6.498 s | 27.282 s |

```text
1T -> 8T speedup: 1.6266x
8T parallel efficiency: 20.33%
```

OpenMP scaling は弱い。

---

# 24. Emilia_923: Stage-4B pipeline

| Threads | Move | Refine | Aggregate | Leiden total | Execution |
|---:|---:|---:|---:|---:|---:|
| 1 | 10.569 | 8.008 | 0.324 | 19.575 | 19.852 |
| 2 | 9.325 | 33.542 | 0.218 | 44.108 | 44.385 |
| 4 | 6.999 | 34.185 | 0.146 | 42.240 | 42.519 |
| 8 | 6.498 | 19.980 | 0.109 | 27.282 | 27.561 |

8T:

```text
RefinePartition ~73.2% of Leiden total
MoveNodesFast   ~23.8%
```

Stage-4B は MoveNodesFast 自体を改善するが、partition trajectory により RefinePartition workload が大きく増える。

---

# 25. Emilia_923: 同一環境 baseline comparison

## Serial

```text
threads        1
MoveNodesFast  4.354 s
Refine         50.184 s
Leiden total   55.504 s
visits         3,742,978
moves          1,307,167
quality        3.59781e7
blocks/colors  37/8
```

## Stage-4A.1

```text
threads        8
MoveNodesFast  12.724 s
Refine          1.108 s
Leiden total   14.617 s
processed      12,377,795
moves           2,563,864
quality         3.60897e7
blocks/colors   33/8
```

## Stage-4B

```text
threads        8
MoveNodesFast   6.498 s
Refine         19.980 s
Leiden total   27.282 s
processed       6,696,856
moves          ~1,341,846
quality median  3.59019e7
blocks          36–40
colors          7–8
```

主要比較:

```text
Stage-4B Move vs Stage-4A.1 Move:
1.958x faster

Stage-4B Move vs serial Move:
0.670x speed
=> Stage-4B is about 1.49x slower than serial

Stage-4B Leiden total vs Stage-4A.1:
0.536x speed
=> Stage-4B pipeline is about 1.87x slower

Stage-4B Leiden total vs serial:
~2.03x faster
```

---

# 26. Emilia_923: Work amplification

| Metric | Stage-4A.1 | Stage-4B 8T |
|---|---:|---:|
| processed/visits | 12,377,795 | 6,696,856 |
| vs serial | 3.307x | 1.789x |
| moves | 2,563,864 | ~1,341,846 |
| move amplification vs serial | 1.961x | 1.027x |

Stage-4B は Stage-4A.1 の work amplification を大幅削減したが、serial より processed vertices は約78.9%多い。

---

# 27. Emilia_923: Stage-4B synchronization profile

8T代表値:

```text
Processed vertices       6,696,856
Successful moves         1,340,601
Failed validations              63
Candidate evaluations   52,527,052
Neighbor scans          276,928,193

Enqueue attempts         45,974,675
Successful enqueues      4,621,282
Duplicate suppressions  41,353,393
Duplicate rate              89.95%

Commit rejection rate       0.00470%

Lock attempts             2,681,328
Lock contentions             11,751
Lock contention rate          0.438%

Empty attempts            6,696,807
Empty claim failures              0

Final rebuild                0.120 s
```

thread-time:

```text
Queue management        31.786 s
Neighbor construction    6.648 s
Candidate evaluation     5.080 s
Affected/requeue update  3.630 s
Commit/revalidation      0.547 s
Lock wait                0.101 s
Final rebuild            0.120 s
```

これらは worker ごとの wall-duration の合計であり、互いに重なっている。wall time と単純加算しない。

---

# 28. Stage-4B 最大の MoveNodesFast bottleneck

## 1位: single shared work queue

```text
8T queue thread-time = 31.786 s
```

## 2位: duplicate activation

```text
enqueue attempts        45.97M
successful enqueue       4.62M
duplicate suppressions  41.35M
duplicate rate          89.95%
```

## 3位: work amplification

```text
processed vs serial = 1.789x
neighbor scans       = 276.9M
candidate evals      = 52.5M
```

一方、以下は主要 bottleneck ではない。

```text
commit success ~99.995%
lock contention ~0.438%
lock wait very small
empty claim failure = 0
```

したがって source/target community lock 最適化は優先度が低い。

---

# 29. Verification sweep

Emilia 8Tでは:

```text
verification sweeps = 20
```

Stage-4B.3 時点では verification 専用 timer がなく、verification overhead fraction は未計測。次段階で独立計測する。

---

# 30. Emilia_923: Quality variability

Stage-4B 8T:

```text
quality median  3.59019e7
range           3.58718e7 – 3.59639e7
mean            3.591042e7
SD              ~30,772

relative vs Stage-4A.1: -0.520%
relative vs serial:     -0.212%

blocks: 36–40
colors: 7–8
levels: 8–9
moves: 1,323,902–1,350,305
```

1T は同じ結果を再現。2/4/8T は async scheduling により run ごとに `.blk` / `.bcol` hash が異なる。これは nondeterministic asynchronous algorithm と整合する。

---

# 31. af_shell3: Stage-4B

| Threads | Move median | Leiden median |
|---:|---:|---:|
| 1 | 2.969 s | 7.578 s |
| 2 | 2.994 s | 5.787 s |
| 4 | 2.232 s | 3.938 s |
| 8 | 2.114 s | 3.223 s |

比較:

```text
Serial Move        1.160 s
Serial Leiden     19.272 s

Stage-4A.1 8T Move   1.488 s
Stage-4A.1 8T Leiden 2.851 s

Stage-4B 8T Move     2.114 s
Stage-4B 8T Leiden   3.223 s
```

Stage-4B は af_shell3 では Stage-4A.1 より遅い。8T profileでも queue thread-time が最大で、lock contention / commit rejection は極小。

---

# 32. Stage-4B.3 最終評価

```text
Stage-4B.3 performance evaluation:
complete

Stage-4B performance:
improved but not yet competitive

OpenMP scaling:
poor

Production default:
Stage-4A.1
```

重要な解釈:

- Stage-4B は Stage-4A.1 の work amplification を削減した。
- Stage-4B の MoveNodesFast は Stage-4A.1 より大幅高速。
- しかし serial MoveNodesFast よりまだ遅い。
- shared queue と duplicate activation が scaling を制限。
- Emilia では Stage-4B が生成する partition が RefinePartition に非常に不利。
- pipeline 全体は Stage-4A.1 より遅い。
- af_shell3 では Refine は問題にならないが、Stage-4B MoveNodesFast overhead が大きい。

---

# 33. 現在の最重要な研究上の知見

```text
並列 local-moving の実行方式は、MoveNodesFast 自身の時間だけでなく、
その後の RefinePartition workload を大きく変える。
```

Stage-4B は MoveNodesFast の work amplification を削減したが、partition trajectory が変わり、Emilia では downstream RefinePartition が増大した。

研究用の表現例:

```text
An asynchronous local-moving scheme reduced the work amplification
of the round-based parallel method, but its different partition
trajectory substantially increased the cost of the subsequent
refinement phase for some graphs.
```

---

# 34. 次の段階: Stage-4B.4

推奨名称:

```text
Stage-4B.4:
Targeted Optimization and Downstream-Work Analysis
```

## 34.1 Queue contention / duplicate activation

最優先。

現状:

```text
single shared mutex queue
~90% duplicate enqueue attempts
```

候補:

- thread-local queues
- per-thread deque
- lightweight work stealing
- activation batching
- duplicate-attempt削減
- neighbor activation policy改善

最初から複雑な custom lock-free queue は導入しない。

## 34.2 Verification sweep 独立計測

追加すべき:

```text
verification sweep count
verification wall time
processed vertices per sweep
successful moves per sweep
candidate evaluations per sweep
```

20 sweeps がどれだけ work amplification に寄与しているか確認する。

## 34.3 Stage-4B partition -> RefinePartition workload 分析

Emilia では非常に重要。

Stage-4A.1 と Stage-4B で level ごとに比較:

```text
graph vertices
graph edges
parent community count
subset count
subset size distribution
max subset size
candidate vertices
processed candidates
refinement moves
sparse target count
delta evaluations
largest subset time
RefinePartition time
thread load balance
```

目的:

```text
なぜ Stage-4A.1 Refine ~1.1 s に対し
Stage-4B Refine ~20 s になるのか
```

を構造的に説明する。

## 34.4 Small/coarse-level serial fallback

coarse level では parallel overhead が支配する可能性がある。

候補:

```text
if n < threshold:
    use serial MoveNodesFast
else:
    Stage-4B
```

threshold は先に決めず、level profiling 後に判断する。

---

# 35. Stage-4B.4で優先すべきでないもの

現時点の profile から、以下は優先度が低い。

- community lock micro-optimization
- source/target lock redesign
- empty community claim optimization
- commit validation optimization

理由:

```text
commit success ~99.995%
lock contention < 0.5%
lock wait small
empty claim failure = 0
```

---

# 36. Profiling 上の注意

Stage-4B.3 時点で `profile macro OFF` でも worker 内の per-thread counter / timer 取得が常時走っている。

したがって primary wall time に instrumentation overhead が含まれている可能性がある。

次の性能最適化段階では `true zero-overhead profile OFF path` を作って比較する価値がある。ただし algorithm 変更と profiling removal を同時に行い、効果を混同しないこと。

---

# 37. 現在の production 方針

```text
Production default:
Stage-4A.1
```

Stage-4B は experimental。

理由:

- MoveNodesFast は改善したが serial より遅い。
- OpenMP scaling が弱い。
- Emilia pipeline が Stage-4A.1 より大幅に遅い。
- quality がわずかに低下・変動。
- queue contention が未解決。

Stage-4B を default に変更してはならない。

---

# 38. Correctness / sanitizer 状況まとめ

```text
normal tests:
pass

Stage-4B dedicated tests:
pass

ASan:
pass

UBSan:
pass

std::thread TSan:
clean

OpenMP/libgomp TSan:
region/lifetime warnings

minimal safe OpenMP program:
equivalent warnings reproduced

OpenMP stress:
1/2/4/8T, 1000 suites, pass

queue/outstanding invariants:
verified

empty-community concurrency:
verified

termination:
verified
```

GCC/libgomp TSan warning は `likely GCC/libgomp–TSan interoperability issue` と扱う。「confirmed false positive」と断定しない。

---

# 39. Cursor / GPT-5.6 Sol への引き継ぎ指示

新しい Cursor セッションでは、まずこの Markdown を読み、以下を守ること。

1. すでに完了済みの逐次最適化を再提案しない。
2. `num_edges()` self-loop 問題は修正済み。
3. AggregateGraph Stage 2 は完成済み。
4. RefinePartition の exact sparse target は完成済み。
5. incremental active-list は production path で使用済み。
6. Stage-4A.1 は production default。
7. Stage-4B は correctness validation 済みだが experimental。
8. Stage-4B の主要 Move bottleneck は queue + duplicate activation。
9. community lock / commit conflict は主要因ではない。
10. Emilia では downstream RefinePartition workload 増大が極めて重要。
11. 次に単純な queue micro-optimizationだけをしない。
12. Stage-4B.4では queue と downstream Refine structure を両方分析する。
13. performance変更前後で correctness / quality を必ず確認する。
14. asynchronous Stage-4B は byte-identical reproducibility を要求しない。
15. benchmark は同一 source/build/environment で比較する。
16. profile timer が nested/thread-time の場合、wall time と単純加算しない。

---

# 40. Cursorへの最初の推奨プロンプト

```text
このプロジェクトの引き継ぎ資料
"MyLeiden 並列化・最適化 開発履歴 / Cursor 引き継ぎノート"
を最初に読み、これまでに完了した最適化・検証結果を前提として作業してください。

現在の branch は Parallel4 です。
production default は Stage-4A.1、
Stage-4B は experimental です。

次の作業は Stage-4B.4:
Targeted Optimization and Downstream-Work Analysis とします。

最初にコードを変更せず、現在の Stage-4B 実装と profiling を監査し、

1. shared work queue contention
2. duplicate activation (~90%)
3. verification sweeps
4. Stage-4B partition が Emilia_923 の RefinePartition を
   Stage-4A.1 の約1.1 sから約20 sへ増大させる原因

を、どの計測を追加すれば分離できるか提案してください。

既存の correctness semantics、
QualityFunction abstraction、
incremental active-list、
exact sparse refinement targets、
AggregateGraph Stage 2
は変更しないでください。

まず分析計画だけを示し、
この段階では実装を変更しないでください。
```

---

# 41. 重要な数値の短縮版

## Emilia_923

```text
Serial:
Move 4.354
Refine 50.184
Leiden 55.504
visits 3.743M
moves 1.307M
quality 3.59781e7
37 blocks / 8 colors

Stage-4A.1 8T:
Move 12.724
Refine 1.108
Leiden 14.617
processed 12.378M
moves 2.564M
quality 3.60897e7
33 / 8

Stage-4B 8T:
Move 6.498
Refine 19.980
Leiden 27.282
processed 6.697M
moves ~1.342M
quality median 3.59019e7
blocks 36–40
colors 7–8
```

Stage-4B:

```text
Move vs Stage-4A.1: 1.958x faster
Move vs serial: 1.49x slower
Leiden vs Stage-4A.1: 1.87x slower
Leiden vs serial: ~2.03x faster
```

8T queue:

```text
enqueue attempts        45.97M
successful enqueues      4.62M
duplicates              41.35M
duplicate rate          89.95%
queue thread-time       31.786 s
lock contention          0.438%
commit success          99.995%
```

## af_shell3

```text
Serial Move     1.160
Serial Leiden  19.272

Stage-4A.1 8T:
Move   1.488
Leiden 2.851

Stage-4B 8T:
Move   2.114
Leiden 3.223
```

Stage-4B is slower than Stage-4A.1 on af_shell3.

---

# 42. 現時点の結論

```text
Stage-4A.1:
work amplificationは大きいが、
EmiliaではRefinePartitionに非常に有利なpartitionを生成し、
pipeline全体では現時点の最良parallel variant。

Stage-4B:
asynchronous化によりStage-4A.1のwork amplificationを大幅に削減し、
MoveNodesFast自体は約2倍高速化した。

しかし、
single shared queue + duplicate activationによりscalingが弱く、
serial MoveNodesFastを超えられていない。

さらにEmiliaではpartition trajectoryの違いにより
RefinePartitionが約1.1 s -> 約20 sへ悪化し、
pipeline全体ではStage-4A.1より遅い。

したがって、
次の最適化ではMoveNodesFastだけでなく、
MoveNodesFastが生成するpartitionと
RefinePartition workloadの連鎖を同時に分析する必要がある。
```

---

# 43. 次にやること

正式な次段階:

```text
Stage-4B.4:
Targeted Optimization and Downstream-Work Analysis
```

第一段階ではコードを変更せず、計測設計から始める。

優先:

```text
1. queue contention / duplicate activation
2. verification sweep cost
3. level-by-level MoveNodesFast statistics
4. Stage-4A.1 vs Stage-4B RefinePartition input structure
5. largest refinement subset / candidate workload
6. small/coarse-level serial fallback feasibility
```

その結果を見てから Stage-4B.4 の実装変更を決める。
