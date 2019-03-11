#ifndef ANNEALING_H
#define ANNEALING_H

#include <vector>
#include <random>
#include "Graph.hpp"

using namespace std;

#define BATCH_SIZE 10000
#define MAX_BATCH_ATTEMPTS BATCH_SIZE * 10
#define EPSILON 0.005
#define EPSILONP 0.005
#define PHI 0.95

class Annealing
{
private:
  int n;
  Graph G;
  Graph GC;
  pair<double, double> *cities;
  double normalizer;
  double maxD;
  vector<int> s;
  vector<int> minS;
  double minSCost;
  default_random_engine dre;
  uniform_int_distribution<int> uid;

  double costFunction(vector<int> &);
  void computeNormalizer(vector<int> &);
  void computeGComplete();
  double computeBatch(double);
  void thresholdAccepting(double);
  double initialTemperature(double, double);
  double acceptedPercentage(double);
  double binarySearch(double, double, double);
  void createInitialSolution(vector<int> &);
  pair<vector<int>, double> neighbour(double);

  static double naturalD(double, double, double, double);

public:
  Annealing(int);
  void setRandomEngine(int, int);
  void validEdge(int, int, double);
  void addCity(int, pair<double, double>);
  pair<vector<int>, double> computeSolution(vector<int> &);
};

#endif