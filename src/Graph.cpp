#include "Graph.hpp"

Graph::Graph() {}

Graph::Graph(int n)
{
  this->order = n;
  this->G = new double *[n];
  for (int i = 0; i < n; i++)
  {
    this->G[i] = new double[n];
    for (int j = 0; j < n; j++)
      if (i == j)
        this->G[i][j] = 0;
      else
        this->G[i][j] = -1;
  }
}

void Graph::addEdge(int u, int v, double weight) { this->G[u][v] = weight; }

double Graph::getWeight(int u, int v) { return this->G[u][v]; }