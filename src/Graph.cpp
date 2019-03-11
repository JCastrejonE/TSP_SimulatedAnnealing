#include "Graph.hpp"

Graph::Graph() {}

Graph::Graph(int n)
{
  order = n;
  G = new double *[n];
  for (int i = 0; i < n; i++)
  {
    G[i] = new double[n];
    for (int j = 0; j < n; j++)
      if (i == j)
        G[i][j] = 0;
      else
        G[i][j] = -1;
  }
}

void Graph::addEdge(int u, int v, double weight) { G[u][v] = weight; }

double Graph::getWeight(int u, int v) { return G[u][v]; }