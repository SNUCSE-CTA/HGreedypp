CXX = g++
CXXFLAGS = -std=gnu++0x -O3 -Iinclude
SRC_TRIANGLE = TrianglePeel.cpp
TARGET_TRIANGLE = TrianglePeel

# SRC_TRIANGLE = exp.cpp
# TARGET_TRIANGLE = exp

NUM := 1

all: $(TARGET_TRIANGLE)

$(TARGET_TRIANGLE): $(SRC_TRIANGLE)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(TARGET_TRIANGLE) output/txt/*

.PHONY: all clean

run: all
	for i in dataset/*.txt; do \
		filename=$$(basename $$i); \
		output_directory="output/TrianglePeel"; \
		if [ ! -d "$$output_directory" ]; then \
			mkdir -p "$$output_directory"; \
		fi; \
		echo "$$i"; \
		./$(TARGET_TRIANGLE) < "$$i" > "$$output_directory/$$filename"; \
	done
