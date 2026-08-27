CXX ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I. -Icommon
LDFLAGS ?=
LDLIBS ?=

TARGETS := test_quality test_move_nodes_fast test_refinement test_aggregate test_coarse_partition test_leiden test_block_eval LeidenCP LeidenMD eblock
TEST_QUALITY_OBJS := QualityFunction.o tests/test_quality.o
TEST_MOVE_NODES_FAST_OBJS := QualityFunction.o Leiden.o tests/test_move_nodes_fast.o
TEST_REFINEMENT_OBJS := QualityFunction.o Leiden.o tests/test_refinement.o
TEST_AGGREGATE_OBJS := QualityFunction.o Leiden.o tests/test_aggregate.o
TEST_COARSE_PARTITION_OBJS := QualityFunction.o Leiden.o tests/test_coarse_partition.o
TEST_LEIDEN_OBJS := QualityFunction.o Leiden.o tests/test_leiden.o
TEST_BLOCK_EVAL_OBJS := common/BlockIO.o common/Coloring.o tests/test_block_eval.o
LEIDEN_CP_OBJS := QualityFunction.o Leiden.o common/mm_io.o common/BlockIO.o common/Coloring.o leiden_main_cp.o
LEIDEN_MD_OBJS := QualityFunction.o Leiden.o common/mm_io.o common/BlockIO.o common/Coloring.o leiden_main_md.o
EBLOCK_OBJS := common/mm_io.o common/BlockIO.o common/Coloring.o eblock.o
DEPS := $(sort $(TEST_QUALITY_OBJS:.o=.d) $(TEST_MOVE_NODES_FAST_OBJS:.o=.d) $(TEST_REFINEMENT_OBJS:.o=.d) $(TEST_AGGREGATE_OBJS:.o=.d) $(TEST_COARSE_PARTITION_OBJS:.o=.d) $(TEST_LEIDEN_OBJS:.o=.d) $(TEST_BLOCK_EVAL_OBJS:.o=.d) $(LEIDEN_CP_OBJS:.o=.d) $(LEIDEN_MD_OBJS:.o=.d) $(EBLOCK_OBJS:.o=.d))

.PHONY: all test clean

all: $(TARGETS)

test_quality: $(TEST_QUALITY_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_move_nodes_fast: $(TEST_MOVE_NODES_FAST_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_refinement: $(TEST_REFINEMENT_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_aggregate: $(TEST_AGGREGATE_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_coarse_partition: $(TEST_COARSE_PARTITION_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_leiden: $(TEST_LEIDEN_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_block_eval: $(TEST_BLOCK_EVAL_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

LeidenCP: $(LEIDEN_CP_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

LeidenMD: $(LEIDEN_MD_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

eblock: $(EBLOCK_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

leiden_main_cp.o: leiden_main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DLEIDEN_DRIVER_CPM -MMD -MP -c $< -o $@

leiden_main_md.o: leiden_main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DLEIDEN_DRIVER_MODULARITY -MMD -MP -c $< -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: $(TARGETS)
	./test_quality
	./test_move_nodes_fast
	./test_refinement
	./test_aggregate
	./test_coarse_partition
	./test_leiden
	./test_block_eval

clean:
	rm -f $(TARGETS) *.o *.d common/*.o common/*.d tests/*.o tests/*.d

-include $(DEPS)
