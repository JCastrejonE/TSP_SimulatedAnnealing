#ifndef GRAPH_H
#define GRAPH_H

template <int n>
class Graph
{
private:
  int order = n;
  double G[n][n];

public:
  Graph();
  void addEdge(int u, int v, double weight) { G[u][v] = weight; }
  double getWeight(int u, int v) { return G[u][v]; }
};

template <int n>
Graph<n>::Graph()
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      if (i == j)
        G[i][j] = 0;
      else
        G[i][j] = -1;
  }
}

#endif
