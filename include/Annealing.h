#ifndef ANNEALING_H
#define ANNEALING_H

#include <iostream>
#include <vector>
#include <random>
#include "Graph.h"

#define BATCH_SIZE 1000
#define MAX_BATCH_ATTEMPTS BATCH_SIZE * BATCH_SIZE
#define EPSILON 0.001
#define EPSILONP 0.000001
#define PHI 0.9

using namespace std;

template <int n>
class Annealing
{
private:
  static Graph<n> G;
  static Graph<n> GC;
  static pair<double, double> cities[n];
  static double normalizer;
  static double maxD;

  static vector<int> s;

  static double naturalD(double, double, double, double);

  static double computeBatch(double);
  static double acceptedPercentage(double);
  static double binarySearch(double, double, double);

public:
  // Annealing();
  static default_random_engine dre;
  static uniform_int_distribution<int> uid;

  static void computeNormalizer(vector<int>);
  static void computeGComplete();
  static void validEdge(int, int, double);
  static void addCity(int, pair<double, double>);
  static double costFunction(vector<int>);

  static double initialTemperature(double, double);
  static void thresholdAccepting(double);
  static void createInitialSolution(vector<int>);
  static vector<int> neighbour(vector<int>);
};

template <int n>
Graph<n> Annealing<n>::G;

template <int n>
Graph<n> Annealing<n>::GC;

template <int n>
pair<double, double> Annealing<n>::cities[];

template <int n>
double Annealing<n>::normalizer = 0;

template <int n>
double Annealing<n>::maxD = 0;

template <int n>
vector<int> Annealing<n>::s;

template <int n>
default_random_engine Annealing<n>::dre;

template <int n>
uniform_int_distribution<int> Annealing<n>::uid(1, n);

template <int n>
double Annealing<n>::costFunction(vector<int> S)
{
  double sum = 0;
  for (int i = 0; i < S.size() - 1; i++)
  {
    sum += Annealing<n>::GC.getWeight(S[i], S[i + 1]);
  }
  double cost = sum / Annealing<n>::normalizer;
  return cost;
}

template <int n>
void Annealing<n>::validEdge(int u, int v, double weight)
{
  Annealing<n>::G.addEdge(u, v, weight);
  Annealing<n>::G.addEdge(v, u, weight);
}

template <int n>
void Annealing<n>::computeNormalizer(vector<int> S)
{
  vector<double> L;

  int m = S.size();

  for (int i = 0; i < m; i++)
  {
    for (int j = i + 1; j < m; j++)
    {
      double w_ij = Annealing<n>::G.getWeight(S[i], S[j]);
      if (w_ij > 0)
        L.push_back(w_ij);
    }
  }

  sort(L.rbegin(), L.rend());
  Annealing<n>::maxD = L.front();
  // printf("Maxd(S): %2.9f\n", Annealing<n>::maxD);
  // cout << "size(S): " << S.size() << endl;
  // cout << "size(L): " << L.size() << endl;
  // cout << "size(L'): " << S.size() - 1 << endl;
  double normalizer = 0;
  for (int i = 0; i < S.size() - 1; i++)
    normalizer += L[i];
  // printf("norm(S): %2.9f\n", normalizer);
  Annealing<n>::normalizer = normalizer;
}

template <int n>
double Annealing<n>::naturalD(double latU, double longU, double latV, double longV)
{
  latU = latU * M_PI / 180;
  latV = latV * M_PI / 180;
  longU = longU * M_PI / 180;
  longV = longV * M_PI / 180;

  double A = pow(sin((latV - latU) / 2), 2) +
             cos(latU) * cos(latV) *
                 pow(sin((longV - longU) / 2), 2);
  double R = 6373000.0; // m
  double C = 2.0 * atan2(sqrt(A), sqrt(1 - A));

  return R * C;
}

template <int n>
void Annealing<n>::computeGComplete()
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      double actualW = Annealing<n>::G.getWeight(i, j);
      if (actualW == -1)
      {
        double latU = cities[i].first;
        double latV = cities[j].first;
        double longU = cities[i].second;
        double longV = cities[j].second;
        double ws = Annealing<n>::naturalD(latU, longU, latV, longV) * Annealing<n>::maxD;
        Annealing<n>::GC.addEdge(i, j, ws);
      }
      else
      {
        Annealing<n>::GC.addEdge(i, j, actualW);
      }
    }
}

template <int n>
void Annealing<n>::addCity(int i, pair<double, double> coords)
{
  Annealing<n>::cities[i] = coords;
}

template <int n>
double Annealing<n>::computeBatch(double T)
{
  int c = 0;
  int i = 0;
  double r = 0.0;
  while (c < BATCH_SIZE)
  {
    vector<int> sp = Annealing<n>::neighbour(Annealing<n>::s);
    if (Annealing<n>::costFunction(sp) < Annealing<n>::costFunction(Annealing<n>::s) + T)
    {
      Annealing<n>::s = sp;
      c += 1;
      r += Annealing<n>::costFunction(sp);
    }
    i += 1;
    if (i >= MAX_BATCH_ATTEMPTS)
      break;
  }
  printf("%2.9f\n", Annealing<n>::costFunction(Annealing<n>::s));
  double res = (double)r / (double)BATCH_SIZE;
  return res;
}

template <int n>
void Annealing<n>::thresholdAccepting(double T)
{
  double p = 0;
  while (T > EPSILON)
  {
    double q = numeric_limits<double>::max();
    while (p <= q)
    {
      q = p;
      p = Annealing<n>::computeBatch(T);
    }
    T = PHI * T;
  }
  for(auto i: Annealing<n>::s) {
    printf("%d,", i);
  }
  printf("\n");
}

template <int n>
double Annealing<n>::initialTemperature(double T, double P)
{
  double p = Annealing<n>::acceptedPercentage(T);
  double T1;
  double T2;
  if (abs(P - p) <= EPSILONP)
    return T;
  if (p < P)
  {
    while (p < P)
    {
      T = 2 * T;
      p = Annealing<n>::acceptedPercentage(T);
    }
    T1 = T / 2;
    T2 = T;
  }
  else
  {
    while (p > P)
    {
      T = T / 2;
      p = Annealing<n>::acceptedPercentage(T);
    }
    T1 = T;
    T2 = 2 * T;
  }
  double res = Annealing<n>::binarySearch(T1, T2, P);
  printf("Ti: %2.9f\n", res);
  return res;
}

template <int n>
double Annealing<n>::acceptedPercentage(double T)
{
  int c = 0;
  int tries = 0;
  for (int i = 0; i < BATCH_SIZE; i++)
  {
    vector<int> sp = Annealing<n>::neighbour(Annealing<n>::s);
    if (Annealing<n>::costFunction(sp) <= Annealing<n>::costFunction(Annealing<n>::s) + T)
    {
      c += 1;
      Annealing<n>::s = sp;
    }
    tries += 1;
    if (tries >= MAX_BATCH_ATTEMPTS)
      break;
  }
  double res = (double)c / (double)BATCH_SIZE;
  // printf("Accepted: %d Batchsize: %d P: %2.9f\n", c, BATCH_SIZE, res);
  return res;
}

template <int n>
double Annealing<n>::binarySearch(double T1, double T2, double P)
{
  double Tm = (T1 + T2) / 2;
  // printf("Binary T1: %2.4f T2: %2.4f Tm: %2.4f P: %2.4f\n", T1, T2, Tm, P);
  if (T2 - T1 < EPSILONP)
    return Tm;
  double p = Annealing<n>::acceptedPercentage(Tm);
  if (abs(P - p) < EPSILONP)
    return Tm;
  if (p > P)
    return Annealing<n>::binarySearch(T1, Tm, P);
  else
    return Annealing<n>::binarySearch(Tm, T2, P);
}

template <int n>
void Annealing<n>::createInitialSolution(vector<int> S)
{
  Annealing<n>::s = S;
  shuffle(Annealing<n>::s.begin(), Annealing<n>::s.end(), dre);
}

template <int n>
vector<int> Annealing<n>::neighbour(vector<int> s)
{
  vector<int> sp = s;
  int i = Annealing<n>::uid(dre);
  int j;
  do {
    j = Annealing<n>::uid(dre);
  } while(i == j);
  int tmp = sp[i];
  sp[i] = sp[j];
  sp[j] = tmp;
  return sp;
}

#endif