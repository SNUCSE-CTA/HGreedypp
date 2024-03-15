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
Demo: ./TrianglePeel.exe 10 < toy.txt > toy.log
Charalampos E. Tsourakakis, babis@seas.harvard.edu
*/

const int MAXV = 10000000;
const int MAXE = 100000000;

int MAXDEG = 0;
// we shall assume that all vertices are numbered  from 1 to n

int E, V;
int **eList;
int *degrees;
int NumTriangles = 0;

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
  cout << "Elapsed time to read input: "
       << (clock() - startTime) / CLOCKS_PER_SEC << endl;
}

void ElapsedTime(double startTime, char *s) {
  cout << "Elapsed time for " << s << ": "
       << (clock() - startTime) / CLOCKS_PER_SEC << endl;
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

int triangles[MAXV + 1];

void CountTrianglesNaively() {
  MaximumDegree();
  for (int i = 1; i <= V; ++i) {
    for (int j1 = 0; j1 < adjacent[i].size(); ++j1) {
      for (int j2 = 0; j2 < j1; ++j2) {
        if (HASH::find(adjacent[i][j1], adjacent[i][j2])) {
          triangles[i]++;
          triangles[adjacent[i][j1]]++;
          triangles[adjacent[i][j2]]++;
          NumTriangles++;
        }
      }
    }
  }
  cout << "Total number of triangles " << NumTriangles << endl;
}

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

double TRIANGLEDENSITY = -1;
int TrianglePeelSize = -1;
double TrianglePeelSizeFraction = -1.0;
int l[MAXV + 1];

int triangles_current[MAXV + 1];
int degrees_current[MAXV + 1];

int OPTIND;

priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>>
    q;

void PQPeelTriangles() {
  TRIANGLEDENSITY = NumTriangles / V;
  TrianglePeelSize = V;
  TrianglePeelSizeFraction = 1.0;

  map<vector<int>, int> D;

  cout << "************** Triangle  PEELING ***************\n";
  cout << setw(9) << "Iteration" << setw(15) << "Elapsed time" << setw(15)
       << "Avg Td" << setw(9) << "Size" << setw(15) << "Size/|V|" << setw(25)
       << "Stopped window size" << endl;

  double startTime = clock();
  int it_cnt = 0;

  vector<int> stop_sign = {3, 5, 10, 20};

  while (true) {
    double numedges = (double)E;
    double numoftriangles = NumTriangles;
    double numvertices = (double)V;
    double max_density = 0;

    int max_density_ind = 0;

    bool densest_found = false;
    bool stop_sign_accepted = false;

    while (!q.empty()) q.pop();
    for (int i = 1; i <= V; ++i) degrees_current[i] = degrees[i];
    for (int i = 1; i <= V; ++i) triangles_current[i] = triangles[i];

    for (int i = 1; i <= V; ++i)
      q.push(make_pair(l[i] + triangles_current[i], i));

    int counter = 1;
    permutation.resize(1);
    while (true) {
      int c = q.top().second;
      q.pop();
      if (triangles_current[c] < 0)  // == -1
        continue;
      l[c] += triangles_current[c];

      permutation.push_back(c);
      counter++;

      numedges -=
          degrees_current[c];  // number of edges goes down by degree of c
      --numvertices;           // one vertex less now
      numoftriangles -= triangles_current[c];
      for (int i = 0; i < AdjList[c].size(); i++) {
        int w = AdjList[c][i];
        degrees_current[w]--;
        if (triangles_current[c] > 0) {
          int Tw = triangles_current[w];
          if (Tw > 0 && degrees_current[c] >= 2) {
            for (int j = 0; j < adjacent[w].size(); j++) {
              int candidate = adjacent[w][j];
              if (triangles_current[candidate] >= 0 && candidate != c) {
                if (degrees_current[c] <= degrees_current[candidate]) {
                  for (int rc = 0; rc < AdjList[c].size(); rc++) {
                    if (AdjList[c][rc] == candidate) {
                      triangles_current[w]--;
                      triangles_current[candidate]--;
                      break;
                    }
                  }
                } else {
                  for (int rc = 0; rc < AdjList[candidate].size(); rc++) {
                    if (AdjList[candidate][rc] == c) {
                      triangles_current[w]--;
                      triangles_current[candidate]--;
                      break;
                    }
                  }
                }
              }
            }

            q.push(make_pair(triangles_current[w] + l[w], w));
          }
        }
      }
      if (numvertices > 0) {
        if (max_density <= numoftriangles / numvertices) {
          max_density = numoftriangles / numvertices;
          max_density_ind = counter;
          if (TRIANGLEDENSITY < max_density) {
            TRIANGLEDENSITY = max_density;
            TrianglePeelSize = numvertices;
            TrianglePeelSizeFraction = numvertices / V;
            OPTIND = counter;
            densest_found = true;
          }
        }
      }
      triangles_current[c] = -1;  // no need to look at it again
      if (counter == V + 1) break;
    }

    cout << setw(9) << it_cnt + 1;
    cout << setprecision(6) << fixed << setw(15)
         << (clock() - startTime) / CLOCKS_PER_SEC;
    cout << setprecision(6) << fixed << setw(15) << TRIANGLEDENSITY;
    cout << setprecision(3) << fixed << setw(9) << TrianglePeelSize;
    cout << setprecision(3) << fixed << setw(15) << TrianglePeelSizeFraction;

    if (densest_found) {
      D.clear();
      vector<int> d(permutation.begin() + OPTIND, permutation.end());
      sort(d.begin(), d.end());
      ++D[d];
    } else if (max_density == TRIANGLEDENSITY) {
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
    it_cnt++;
  }
}

int main(int argc, char **argv) {
  double startTime = clock();
  GraphIn();
  ElapsedTime(startTime);
  SimplifyGraph();
  BuildGraph();
  BuildAdjList();

  cout << "**************** Triangle  Counting *****************" << endl;
  startTime = clock();
  CountTrianglesNaively();
  ElapsedTime(startTime, "Triangle Counting");
  Deallocate();

  PQPeelTriangles();

  return 0;
}
