CXX ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I. -Icommon
LDFLAGS ?=
LDLIBS ?=

TARGETS := test_quality test_move_nodes_fast test_refinement test_aggregate test_coarse_partition
TEST_QUALITY_OBJS := QualityFunction.o test_quality.o
TEST_MOVE_NODES_FAST_OBJS := QualityFunction.o Leiden.o test_move_nodes_fast.o
TEST_REFINEMENT_OBJS := QualityFunction.o Leiden.o test_refinement.o
TEST_AGGREGATE_OBJS := QualityFunction.o Leiden.o test_aggregate.o
TEST_COARSE_PARTITION_OBJS := QualityFunction.o Leiden.o test_coarse_partition.o
DEPS := $(sort $(TEST_QUALITY_OBJS:.o=.d) $(TEST_MOVE_NODES_FAST_OBJS:.o=.d) $(TEST_REFINEMENT_OBJS:.o=.d) $(TEST_AGGREGATE_OBJS:.o=.d) $(TEST_COARSE_PARTITION_OBJS:.o=.d))

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

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: $(TARGETS)
	./test_quality
	./test_move_nodes_fast
	./test_refinement
	./test_aggregate
	./test_coarse_partition

clean:
	rm -f $(TARGETS) *.o *.d common/*.o common/*.d

-include $(DEPS)
