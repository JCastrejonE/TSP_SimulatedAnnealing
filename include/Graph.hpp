#ifndef GRAPH_H
#define GRAPH_H

class Graph
{
private:
  int order;
  double **G;

public:
  Graph();
  Graph(int);
  void addEdge(int, int, double);
  double getWeight(int, int);
};

#endif
