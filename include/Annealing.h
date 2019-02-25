#ifndef ANNEALING_H
#define ANNEALING_H

#include <iostream>
#include <vector>
#include "Graph.h"

using namespace std;

template <int n>
class Recocido
{
private:
  static Graph<n> G;
  static Graph<n> GC;
  static pair<double, double> cities[n];
  static double maxD;
  static double calculaNormalizador(vector<int>);
  static void calculaGCompleta();
  static double dNatural(double, double, double, double);

public:
  //Recocido();
  static void aristaValida(int, int, double);
  static void addCity(int, pair<double, double>);
  static double costFunction(vector<int>);
};

template <int n>
Graph<n> Recocido<n>::G;

template <int n>
Graph<n> Recocido<n>::GC;

template <int n>
pair<double, double> Recocido<n>::cities[];

template <int n>
double Recocido<n>::maxD=0;

template <int n>
double Recocido<n>::costFunction(vector<int> S)
{
  double normalizador = Recocido<n>::calculaNormalizador(S);
  Recocido<n>::calculaGCompleta();
  double sum = 0;
  for (int i = 0; i < S.size() - 1; i++)
  {
    sum += Recocido<n>::GC.getWeight(S[i], S[i + 1]);
  }
  double cost = sum / normalizador;
  return cost;
}

template <int n>
void Recocido<n>::aristaValida(int u, int v, double weight)
{
  Recocido<n>::G.addEdge(u, v, weight);
  Recocido<n>::G.addEdge(v, u, weight);
}

template <int n>
double Recocido<n>::calculaNormalizador(vector<int> S)
{
  vector<double> L;

  int m = S.size();

  for (int i = 0; i < m; i++)
  {
    for (int j = i + 1; j < m; j++)
    {
      double w_ij = Recocido<n>::G.getWeight(S[i], S[j]);
      if (w_ij > 0)
        L.push_back(w_ij);
    }
  }

  sort(L.rbegin(), L.rend());
  Recocido<n>::maxD = L.front();
  printf("Maxd(S): %2.9f\n", Recocido<n>::maxD);
  cout << "size(S): " << S.size() << endl;
  cout << "size(L): " << L.size() << endl;
  cout << "size(L'): " << S.size()-1 << endl;
  double normalizador = 0;
  for (int i=0; i<S.size()-1; i++)
    normalizador += L[i];
  printf("norm(S): %2.9f\n", normalizador);
  return normalizador;
}

template <int n>
double Recocido<n>::dNatural(double latU, double longU, double latV, double longV)
{
  latU = latU * M_PI / 180;
  latV = latV * M_PI / 180;
  longU = longU * M_PI / 180;
  longV = longV * M_PI / 180;

  double A = pow(sin((latV - latU) / 2), 2) +
             cos(latU) * cos(latV) *
                 pow(sin((longV - longU) / 2), 2);
  double R = 6373000.0; // Km
  double C = 2.0 * atan2(sqrt(A), sqrt(1 - A));

  return R * C;
}

template <int n>
void Recocido<n>::calculaGCompleta()
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      double actualW = Recocido<n>::G.getWeight(i, j);
      if (actualW == -1)
      {
        double latU = cities[i].first;
        double latV = cities[j].first;
        double longU = cities[i].second;
        double longV = cities[j].second;
        double ws = Recocido<n>::dNatural(latU, longU, latV, longV) * Recocido<n>::maxD;
        Recocido<n>::GC.addEdge(i, j, ws);
      }
      else
      {
        Recocido<n>::GC.addEdge(i, j, actualW);
      }
    }
}

template <int n>
void Recocido<n>::addCity(int i, pair<double, double> coords)
{
  Recocido<n>::cities[i] = coords;
}

#endif