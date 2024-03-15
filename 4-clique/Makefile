CXX = g++
CXXFLAGS = -std=gnu++0x -O3 -Iinclude
SRC_TRIANGLE = TrianglePeel.cpp
TARGET_TRIANGLE = TrianglePeel

NUM := 1

all: $(TARGET_TRIANGLE)

$(TARGET_TRIANGLE): $(SRC_TRIANGLE)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(TARGET_TRIANGLE) output/txt/*

.PHONY: all clean

run: all
	for i in ../Improved/dataset/*.txt; do \
		filename=$$(basename $$i); \
		output_directory="../Improved/output/4-clique/"; \
		if [ ! -d "$$output_directory" ]; then \
			mkdir -p "$$output_directory"; \
		fi; \
		echo "$$i"; \
		./$(TARGET_TRIANGLE) $(NUM) < "$$i" > "$$output_directory/$$filename"; \
	done
