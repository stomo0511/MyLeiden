CXX ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I. -Icommon
LDFLAGS ?=
LDLIBS ?=

TARGETS := test_quality
TEST_QUALITY_OBJS := QualityFunction.o test_quality.o
DEPS := $(TEST_QUALITY_OBJS:.o=.d)

.PHONY: all test clean

all: test_quality

test_quality: $(TEST_QUALITY_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: test_quality
	./test_quality

clean:
	rm -f $(TARGETS) *.o *.d common/*.o common/*.d

-include $(DEPS)
