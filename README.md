# Experiment Methods

## 1. Running Experiments with `3-clique`
Navigate to the `3-clique` directory and execute:
```
make run
```
This command processes six datasets located in `3-clique/dataset/` and saves the results in `3-clique/output/TrianglePeel/`. The results show when the process stops for each window size (3, 5, 10, 20). To modify the stopping condition, adjust the `stop_sign` value:

- For `3-clique/TrianglePeel.cpp`, edit line 184.
- For `4-clique/TrianglePeel.cpp`, edit line 187.

## 2. Running Experiments with `4-clique`
Navigate to the `4-clique` directory and execute:
```
make run
```
This also processes the six datasets in `3-clique/dataset/` and saves the results in `3-clique/output/4-clique/`. Similar to `3-clique`, the stopping points for window sizes (3, 5, 10, 20) are recorded.
