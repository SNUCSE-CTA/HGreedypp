#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#include "IO.h"
#include "hash.h"

#define MAX(a, b) (((a) >= (b)) ? (a) : (b))

using namespace std;

/*
Compile: g++ TrianglePeel.cpp -o TrianglePeel -std=gnu++0x -O3
Demo: ./TrianglePeel.exe < toy.txt > toy.log
Charalampos E. Tsourakakis, babis@seas.harvard.edu
*/

const int MAXV = 10000000;
const int MAXE = 100000000;
const int QUERYBUFFER = 6000000;

int MAXDEG = 0;
// we shall assume that all vertices are numbered  from 1 to n

int E, V;
int **eList;
int *degrees;
int NumTriangles = 0;
int NumCliques = 0;

vector<vector<int>> AdjList;
vector<vector<int>> adjacent;

void SimplifyGraph() {
  int E1 = 1;
  int multipleEdges = 0;
  HASH::init();
  for (int i = 1; i <= E; ++i) {
    if (1 <= eList[0][i] && eList[0][i] <= V && 1 <= eList[1][i] &&
        eList[1][i] <= V && eList[0][i] != eList[1][i] &&
        HASH::insert(eList[0][i], eList[1][i])) {
      eList[0][E1] = eList[0][i];
      eList[1][E1] = eList[1][i];
      degrees[eList[0][i]]++;
      degrees[eList[1][i]]++;
      E1++;
    } else
      multipleEdges++;
  }
  E = E1 - 1;

  cout << "Number of edges in simple graph:" << E << endl;
  cout << "Number of multiple edges and self-loops that were removed:"
       << multipleEdges << endl;
}

void GraphIn() {
  int u, v;
  cin >> V >> E;
  cout << "Number of vertices and edges (" << V << "," << E << ")" << endl;
  degrees = new int[V + 1];
  eList = new int *[2];
  for (int i = 0; i < 2; i++) {
    eList[i] = new int[E + 1];
  }
  for (int i = 1; i <= E; ++i) {
    cin >> u >> v;
    if (v > u) {
      eList[0][i] = u;
      eList[1][i] = v;
    }
    if (u > v) {
      eList[0][i] = v;
      eList[1][i] = u;
    }
  }
  AdjList.resize(V + 1);
  adjacent.resize(V + 1);
}

void MaximumDegree() {
  for (int i = 1; i <= V; i++)
    if (MAXDEG < degrees[i]) MAXDEG = degrees[i];
}

void ElapsedTime(double startTime) {
  printf("Elapsed time to read input: %f\n",
         (clock() - startTime) / CLOCKS_PER_SEC);
}

void ElapsedTime(double startTime, char *s) {
  printf("Elapsed time for %s: %f\n", s,
         (clock() - startTime) / CLOCKS_PER_SEC);
}
///////////////adjacency list
int *eStart;
int *mynext;

const int NOEDGE = -1;

void BuildGraph() {
  eStart = new int[V + 1];
  mynext = new int[E + 1];
  for (int i = 1; i <= V; ++i) {
    eStart[i] = NOEDGE;
  }
  for (int i = 1; i <= E; ++i) {
    mynext[i] = eStart[eList[0][i]];
    eStart[eList[0][i]] = i;
  }
}

int cliques[MAXV + 1];

void CountTrianglesNaively() {
  MaximumDegree();
  for (int i = 1; i <= V; ++i) {
    for (int j1 = 0; j1 < adjacent[i].size(); ++j1) {
      for (int j2 = 0; j2 < j1; ++j2) {
        if (HASH::find(adjacent[i][j1], adjacent[i][j2])) {
          NumTriangles += 1;
          for (int j3 = 0; j3 < j2; ++j3) {
            if (HASH::find(adjacent[i][j1], adjacent[i][j3]) &&
                HASH::find(adjacent[i][j2], adjacent[i][j3])) {
              cliques[i]++;
              cliques[adjacent[i][j1]]++;
              cliques[adjacent[i][j2]]++;
              cliques[adjacent[i][j3]]++;
              NumCliques++;
            }
          }
        }
      }
    }
  }
  printf("Total number of triangles %d\n", NumTriangles);
  printf("Total number of 4-cliques %d\n", NumCliques);
}

int MAXTRI = -1;

void BuildAdjList() {
  for (int i = 1; i <= E; i++) {
    AdjList[eList[0][i]].push_back(eList[1][i]);
    AdjList[eList[1][i]].push_back(eList[0][i]);
  }
  for (int i = 1; i <= V; ++i) {
    for (int j = eStart[i]; j != NOEDGE; j = mynext[j]) {
      adjacent[i].push_back(eList[1][j]);
    }
  }
}

void Deallocate() {
  for (int i = 0; i < 2; i++) {
    delete[] eList[i];
  }
  delete[] eList;
  delete[] eStart;
  delete[] mynext;
}

vector<int> permutation;

double CLIQUEDENSITY = -1;
int CliquePeelSize = -1;
double CliqueOptVals[MAXV + 1];
int l[MAXV + 1];
int cliques_current[MAXV + 1];
int degrees_current[MAXV + 1];

int OPTIND;

priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>>
    q;

void PQPeelTriangles() {
  CLIQUEDENSITY = NumCliques / V;
  CliquePeelSize = V;

  map<vector<int>, int> D;

  double startTime = clock();

  vector<int> stop_sign = {3, 5, 10, 20};

  cout << "************** 4-Clique  PEELING ***************" << endl;
  cout << setw(9) << "Iteration" << setw(15) << "Elapsed time" << setw(15)
       << "Avg Cd" << setw(9) << "Size" << setw(15) << "Size/|V|" << setw(25)
       << "Stopped window size" << endl;

  for (int it_cnt = 0; it_cnt < 10000000; it_cnt++) {
    double numedges = (double)E;
    double numofcliques = NumCliques;
    double numvertices = (double)V;

    double max_density = 0;
    int max_density_ind = 0;
    bool densest_found = false;
    bool stop_sign_accepted = false;

    while (!q.empty()) q.pop();
    for (int i = 1; i <= V; ++i) degrees_current[i] = degrees[i];
    for (int i = 1; i <= V; ++i) cliques_current[i] = cliques[i];

    for (int i = 1; i <= V; ++i)
      q.push(make_pair(l[i] + cliques_current[i], i));

    int counter = 1;
    permutation.resize(1);
    while (!q.empty()) {
      int c = q.top().second;
      q.pop();
      if (cliques_current[c] < 0)  // == -1
        continue;
      l[c] += cliques_current[c];

      CliqueOptVals[counter] = numofcliques / numvertices;
      permutation.push_back(c);

      numedges -=
          degrees_current[c];  // number of edges goes down by degree of c
      --numvertices;           // one vertex less now
      numofcliques -= cliques_current[c];
      for (int i = 0; i < AdjList[c].size(); i++) {
        int w = AdjList[c][i];
        degrees_current[w]--;
        if (cliques_current[c] > 0) {
          int Tw = cliques_current[w];
          if (Tw >= 0 && degrees_current[c] >= 3) {
            for (int j = 0; j < adjacent[w].size(); j++) {
              int candidate = adjacent[w][j];
              if (cliques_current[candidate] >= 0 && candidate != c &&
                  HASH::find(candidate, c)) {
                for (int k = 0; k < adjacent[candidate].size(); k++) {
                  int cand = adjacent[candidate][k];
                  if (cliques_current[cand] >= 0 && cand != c && cand != w &&
                      HASH::find(cand, c) && HASH::find(cand, w)) {
                    cliques_current[w]--;
                    cliques_current[candidate]--;
                    cliques_current[cand]--;
                  }
                }
              }
            }

            q.push(make_pair(cliques_current[w] + l[w], w));
          }
        }
      }
      cliques_current[c] = -1;  // no need to look at it again
      if (numvertices > 0) {
        if (max_density < numofcliques / numvertices) {
          max_density = numofcliques / numvertices;
          max_density_ind = counter;
          if (CLIQUEDENSITY < max_density) {
            CLIQUEDENSITY = max_density;
            CliquePeelSize = numvertices;
            OPTIND = counter;
            densest_found = true;
          }
        }
      }
      if (counter == V) break;
    }
    cout << setw(9) << it_cnt + 1;
    printf("%15f", (clock() - startTime) / CLOCKS_PER_SEC);
    cout << setprecision(6) << fixed << setw(15) << CLIQUEDENSITY;
    cout << setprecision(3) << fixed << setw(9) << CliquePeelSize;
    cout << setprecision(3) << fixed << setw(15) << 0;

    if (densest_found) {
      D.clear();
      vector<int> d(permutation.begin() + OPTIND, permutation.end());
      sort(d.begin(), d.end());
      ++D[d];
    } else if (max_density == CLIQUEDENSITY) {
      vector<int> d(permutation.begin() + OPTIND, permutation.end());
      sort(d.begin(), d.end());
      if (++D[d] == stop_sign[0]) {
        cout << setprecision(3) << fixed << setw(25) << stop_sign[0] << endl;
        stop_sign.erase(stop_sign.begin());
        if (stop_sign.size() == 0) break;
        stop_sign_accepted = true;
      }
    }
    if (!stop_sign_accepted) cout << endl;
  }
}

int main(int argc, char **argv) {
  double startTime = clock();
  GraphIn();
  ElapsedTime(startTime);
  SimplifyGraph();
  BuildGraph();
  BuildAdjList();

  /* Count first */
  cout << "**************** Triangle  Counting *****************" << endl;
  startTime = clock();
  CountTrianglesNaively();
  ElapsedTime(startTime, "Triangle Counting");
  Deallocate();

  PQPeelTriangles();

  return 0;
}
