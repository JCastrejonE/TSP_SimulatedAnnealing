#include <iostream>
#include <vector>
#include <random>
#include "Annealing.hpp"
#include "Graph.hpp"

using namespace std;

Annealing::Annealing(int size)
{
  this->n = size;
  this->G = Graph(size);
  this->GC = Graph(size);
  this->cities = new pair<double, double>[size];
}

pair<vector<int>, double> Annealing::computeSolution(vector<int> &S, bool sweep)
{
  this->normalizer = 0;
  this->maxD = 0;
  this->minSCost = numeric_limits<double>::max();

  this->computeNormalizer(S);
  this->computeGComplete();
  this->createInitialSolution(S);
  double Ti = this->initialTemperature(8, .95);
  this->thresholdAccepting(Ti);
  if (sweep)
  {
    this->computeSweep();
  }
  return {minS, minSCost};
}

void Annealing::computeSweep()
{
  double minCost = this->minSCost * this->normalizer;
  double actualCost = minCost;
  bool improved = true;

  while (improved)
  {
    improved = false;
    this->s = this->minS;
    for (int i = 0; i < this->s.size(); i++)
    {
      for (int j = i + 1; j < this->s.size(); j++)
      {
        if (i != j - 1)
        {
          actualCost -= this->GC.getWeight(this->s[i], this->s[i + 1]);
          actualCost -= this->GC.getWeight(this->s[j - 1], this->s[j]);
          actualCost += this->GC.getWeight(this->s[j], this->s[i + 1]);
          actualCost += this->GC.getWeight(this->s[j - 1], this->s[i]);
        }
        if (i != 0)
        {
          actualCost -= this->GC.getWeight(this->s[i - 1], this->s[i]);
          actualCost += this->GC.getWeight(this->s[i - 1], this->s[j]);
        }
        if (j != this->s.size() - 1)
        {
          actualCost -= this->GC.getWeight(this->s[j], this->s[j + 1]);
          actualCost += this->GC.getWeight(this->s[i], this->s[j + 1]);
        }
        this->s[i] = this->s[i] + this->s[j];
        this->s[j] = this->s[i] - this->s[j];
        this->s[i] = this->s[i] - this->s[j];
        if(actualCost < minCost) {
          improved = true;
          minCost = actualCost;
          this->minS = this->s;
        }
      }
    }
  }
  this->minSCost = minCost / this->normalizer;
}

double Annealing::costFunction(vector<int> &S)
{
  double sum = 0;
  for (int i = 0; i < S.size() - 1; i++)
  {
    sum += this->GC.getWeight(S[i], S[i + 1]);
  }
  double cost = sum / this->normalizer;
  return cost;
}

void Annealing::validEdge(int u, int v, double weight)
{
  this->G.addEdge(u, v, weight);
  this->G.addEdge(v, u, weight);
}

void Annealing::computeNormalizer(vector<int> &S)
{
  vector<double> L;

  int m = S.size();

  for (int i = 0; i < m; i++)
  {
    for (int j = i + 1; j < m; j++)
    {
      double w_ij = this->G.getWeight(S[i], S[j]);
      if (w_ij > 0)
        L.push_back(w_ij);
    }
  }

  sort(L.rbegin(), L.rend());
  this->maxD = L.front();
  double normalizer = 0;
  for (int i = 0; i < S.size() - 1; i++)
    normalizer += L[i];
  this->normalizer = normalizer;
}

double Annealing::naturalD(double latU, double longU, double latV, double longV)
{
  latU = latU * M_PI / 180;
  latV = latV * M_PI / 180;
  longU = longU * M_PI / 180;
  longV = longV * M_PI / 180;

  double A = pow(sin((latV - latU) / 2), 2) +
             cos(latU) * cos(latV) *
                 pow(sin((longV - longU) / 2), 2);
  double R = 6373000.0; // meters
  double C = 2.0 * atan2(sqrt(A), sqrt(1 - A));

  return R * C;
}

void Annealing::computeGComplete()
{
  for (int i = 0; i < this->n; i++)
    for (int j = 0; j < this->n; j++)
    {
      double actualW = this->G.getWeight(i, j);
      if (actualW == -1)
      {
        double latU = cities[i].first;
        double latV = cities[j].first;
        double longU = cities[i].second;
        double longV = cities[j].second;
        double ws = Annealing::naturalD(latU, longU, latV, longV) * this->maxD;
        this->GC.addEdge(i, j, ws);
      }
      else
      {
        this->GC.addEdge(i, j, actualW);
      }
    }
}

void Annealing::addCity(int i, pair<double, double> coords)
{
  this->cities[i] = coords;
}

double Annealing::computeBatch(double T)
{
  int c = 0;
  int i = 0;
  double r = 0.0;
  double sCost = this->costFunction(this->s);
  while (c < BATCH_SIZE)
  {
    pair<vector<int>, double> neighbour = this->neighbour(sCost);
    vector<int> sp;
    sp.swap(neighbour.first);
    double spCost = neighbour.second;
    if (spCost < sCost + T)
    {
      if (spCost < minSCost)
      {
        this->minSCost = spCost;
        this->minS = sp;
      }
      sp.swap(this->s);
      sCost = spCost;
      c += 1;
      r += spCost;
      i = 0;
    }
    else
      i += 1;
    if (i >= MAX_BATCH_ATTEMPTS)
      break;
  }
  double res = (double)r / (double)BATCH_SIZE;
  return res;
}

void Annealing::thresholdAccepting(double T)
{
  double p = 0;
  while (T > EPSILON)
  {
    double q = numeric_limits<double>::max();
    while (p <= q)
    {
      q = p;
      p = this->computeBatch(T);
    }
    T = PHI * T;
  }
}

double Annealing::initialTemperature(double T, double P)
{
  double p = this->acceptedPercentage(T);
  double T1;
  double T2;
  if (abs(P - p) <= EPSILONP)
    return T;
  if (p < P)
  {
    while (p < P)
    {
      T = 2 * T;
      p = this->acceptedPercentage(T);
    }
    T1 = T / 2;
    T2 = T;
  }
  else
  {
    while (p > P)
    {
      T = T / 2;
      p = this->acceptedPercentage(T);
    }
    T1 = T;
    T2 = 2 * T;
  }
  double res = this->binarySearch(T1, T2, P);
  return res;
}

double Annealing::acceptedPercentage(double T)
{
  int c = 0;
  int tries = 0;
  for (int i = 0; i < BATCH_SIZE; i++)
  {
    double sCost = this->costFunction(this->s);
    pair<vector<int>, double> neighbour = this->neighbour(sCost);
    vector<int> sp;
    sp.swap(neighbour.first);
    double spCost = neighbour.second;
    if (spCost <= sCost + T)
    {
      c += 1;
      sp.swap(this->s);
      tries = 0;
    }
    else
      tries += 1;
    if (tries >= MAX_BATCH_ATTEMPTS)
      break;
  }
  double res = (double)c / (double)BATCH_SIZE;
  return res;
}

double Annealing::binarySearch(double T1, double T2, double P)
{
  double Tm = (T1 + T2) / 2;
  if (T2 - T1 < EPSILONP)
    return Tm;
  double p = this->acceptedPercentage(Tm);
  if (abs(P - p) < EPSILONP)
    return Tm;
  if (p > P)
    return this->binarySearch(T1, Tm, P);
  else
    return this->binarySearch(Tm, T2, P);
}

void Annealing::createInitialSolution(vector<int> &S)
{
  this->s = S;
  shuffle(this->s.begin(), this->s.end(), dre);
}

pair<vector<int>, double> Annealing::neighbour(double sCost)
{
  vector<int> sp = this->s;
  double spCost = sCost * this->normalizer;
  int i = this->uid(dre);
  int j;
  do
  {
    j = this->uid(dre);
  } while (i == j);

  if (j < i)
  {
    i = i + j;
    j = i - j;
    i = i - j;
  }

  if (i != j - 1)
  {
    spCost -= this->GC.getWeight(sp[i], sp[i + 1]);
    spCost -= this->GC.getWeight(sp[j - 1], sp[j]);
    spCost += this->GC.getWeight(sp[j], sp[i + 1]);
    spCost += this->GC.getWeight(sp[j - 1], sp[i]);
  }

  if (i != 0)
  {
    spCost -= this->GC.getWeight(sp[i - 1], sp[i]);
    spCost += this->GC.getWeight(sp[i - 1], sp[j]);
  }
  if (j != sp.size() - 1)
  {
    spCost -= this->GC.getWeight(sp[j], sp[j + 1]);
    spCost += this->GC.getWeight(sp[i], sp[j + 1]);
  }

  sp[i] = sp[i] + sp[j];
  sp[j] = sp[i] - sp[j];
  sp[i] = sp[i] - sp[j];
  return {sp, spCost / this->normalizer};
}

void Annealing::setRandomEngine(int seed, int size)
{
  default_random_engine dre(seed);
  this->dre = dre;
  uniform_int_distribution<int> uid(0, size - 1);
  this->uid = uid;
}