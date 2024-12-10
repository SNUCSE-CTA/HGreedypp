# H-Greedy++

## Experiment Methods

### 1. Running Experiments with `3-clique`
Navigate to the `3-clique` directory and execute:
```
make run
```
This command processes six datasets located in `3-clique/dataset/` and saves the results in `3-clique/output/TrianglePeel/`. The results show when the process stops for each window size (3, 5, 10, 20). To modify the stopping condition, adjust the `stop_sign` value:

- For `3-clique/TrianglePeel.cpp`, edit line 184.
- For `4-clique/TrianglePeel.cpp`, edit line 187.

### 2. Running Experiments with `4-clique`
Navigate to the `4-clique` directory and execute:
```
make run
```
This also processes the six datasets in `3-clique/dataset/` and saves the results in `3-clique/output/4-clique/`. Similar to `3-clique`, the stopping points for window sizes (3, 5, 10, 20) are recorded.

## Reference
[1] E. Fratkin, B. T. Naughton, D. L. Brutlag, and S. Batzoglou. Motifcut: regulatory motifs finding with maximum density subgraphs. Bioinformatics, 22(14):e150–e157, 2006.  
[2] X. Du, R. Jin, L. Ding, V. E. Lee, and J. H. T. Jr. Migration motif: a spatial - temporal pattern mining  approach for financial markets. SIGKDD, pages 1135–1144, 2009.  
[3] A. Angel, N. Sarkas, N. Koudas, and D. Srivastava. Dense subgraph maintenance under streaming edge weight updates for real-time story identification. PVLDB, 5(6):574–585, 2012.  
[4] C.E. Tsourakakis, The K-clique Densest Subgraph Problem. WWW, pages 1122-1132, 2015.  
[5] M. Danisch, T. H. Chan, and M. Sozio. Large scale density-friendly graph decomposition via convex programming. WWW, pages 233–242, 2017.  
[6] Bintao Sun, Maximilien Danisch, T-H. Hubert Chan, and Mauro Sozio. KClist++: a simple algorithm for finding k-clique densest subgraphs in large graphs. PVLDB, 13(10):1628-1640, 2020.  
[7] Digvijay Boob, Yu Gao, Richard Peng, Saurabh Sawlani, Charalampos Tsourakakis, Di Wang, and Junxing Wang. Flowless: Extracting Densest Subgraphs Without Flow Computations. WWW, pages 573–583, 2020.  
[8] Chekuri, Chandra & Quanrud, Kent & Torres, Manuel. Densest Subgraph: Supermodularity, Iterative Peeling, and Flow. SODA, pages 1531-1555, 2022.  
