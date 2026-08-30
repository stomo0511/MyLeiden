# Stage-4B Remaining Issues

This note summarizes the remaining Stage-4B issues after the Stage-4B.4
measurement pass.

## Current Status

- Production default remains Stage-4A.1.
- Stage-4B remains experimental.
- Stage-4B correctness tests pass, but performance is not yet competitive
  enough to become the default.
- Batched enqueue is a promising Stage-4B experiment, guarded by
  `ENABLE_MOVENODESFAST_STAGE4B_BATCHED_ENQUEUE`.

## Recent Batched Enqueue Median Result

Environment: local macOS build, `-O2`, 8 OpenMP threads, profiling enabled.

```text
Emilia_923
  MoveNodesFast : 5.63339 -> 5.40536 s  (-4.05%)
  Refine        : 4.28206 -> 4.11784 s  (-3.84%)
  Leiden total  : 10.2627 -> 9.89573 s  (-3.58%)
  Queue time    : 30.9019 -> 28.9725 s  (-6.24%)
  Quality median: 3.59383e7 -> 3.59545e7

af_shell3
  MoveNodesFast : 1.80921 -> 1.77449 s  (-1.92%)
  Refine        : 2.19202 -> 2.11192 s  (-3.65%)
  Leiden total  : 4.20232 -> 4.03888 s  (-3.89%)
  Queue time    : 10.4564 -> 10.2856 s  (-1.63%)
  Quality median: 1.63083e7 -> 1.63088e7
```

Batched enqueue reduces successful enqueue publication cost, but it does not
materially reduce duplicate activation attempts.

## Remaining Issues

### 1. Duplicate Activation

Stage-4B still spends substantial work attempting to enqueue vertices that are
already queued or processing.

- Duplicate rate remains around 90% on the large first level.
- Batched enqueue does not reduce duplicate attempts.
- The next optimization should target activation policy or requeue policy, not
  only queue publication mechanics.

### 2. Shared Work Queue Contention

The single shared queue remains a major bottleneck.

- Queue thread-time is still the largest Stage-4B profile component.
- Batched enqueue improves queue time modestly.
- Possible next experiments:
  - per-thread queues
  - local queue plus work stealing
  - batch pop and batch push
  - activation batching with duplicate filtering

Avoid starting with a complex custom lock-free queue.

### 3. Downstream RefinePartition Workload

Stage-4B changes the partition trajectory, which changes RefinePartition input
structure.

Important metrics to keep comparing level by level:

- parent community count
- subset size distribution
- maximum subset size
- sparse target iterations
- active community iterations
- delta evaluations
- refinement moves
- thread load balance

MoveNodesFast speed alone is not enough; downstream RefinePartition cost can
erase the gain.

### 4. Verification Sweep Cost

Final verification sweep cost is measurable.

- Emilia_923 median final verification wall time: about 0.54 s
- af_shell3 median final verification wall time: about 0.28 s

This is not the top bottleneck, but it is large enough to track. Future work
should consider whether all full sweeps are necessary, or whether a smaller
verification frontier is possible without weakening correctness.

### 5. Serial MoveNodesFast Gap

Stage-4B still has overhead relative to the serial local-moving path.

Likely contributors:

- atomic label loads
- community mass locks during candidate evaluation
- shared queue management
- repeated full sweep initialization
- duplicate activation work

### 6. Coarse-Level Parallel Overhead

At small/coarse levels, Stage-4B parallel overhead can dominate useful work.

Potential experiment:

```text
if n < threshold:
    use serial MoveNodesFast
else:
    use Stage-4B
```

Do not choose the threshold upfront. Derive it from level-by-level timing and
workload profiles.

### 7. Profile Overhead Separation

Performance comparisons need a true profile-off path.

- Profile-enabled builds are useful for diagnosis.
- Final timing should compare profile-off baseline vs profile-off optimized
  variants.
- Do not mix profiling removal with algorithmic changes when attributing speedup.

### 8. Quality and Reproducibility

Stage-4B is asynchronous and does not promise byte-identical output across
thread schedules.

Required evaluation metrics:

- quality median and range
- block count range
- color count range
- Leiden level count range
- downstream ICCG performance, when available

Median runtime alone is insufficient for production adoption.

## Near-Term Recommendation

Keep batched enqueue as an experimental macro. The next Stage-4B.4 step should
target duplicate activation and queue structure together, preferably with a
simple per-thread queue or batch-pop/batch-push experiment, while continuing to
record RefinePartition detailed profiles.
