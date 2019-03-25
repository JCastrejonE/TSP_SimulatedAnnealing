#include <iostream>
#include <sstream>
#include <chrono>
#include <fstream>
#include <sqlite3.h>
#include "jarro2783/cxxopts.hpp"
#include "Graph.hpp"
#include "Annealing.hpp"

using namespace std;

// DB FILE PATH
#define DBPATH "../db/tsp.db"

// NUMBER OF CITIES TO READ FROM DB
#define N 1092

// GETS EXECUTED FOR EACH ROW OF THE DB QUERY RESULT.
static int callback(void *, int, char **, char **);

int main(int argc, char** argv)
{
  auto gstart = chrono::steady_clock::now();
  int sseed=0, eseed=10;
  istream *is = &cin;
  ifstream inFile;
  bool hybrid=false, sweep=true;

  cxxopts::Options options(argv[0], "TSP problem solution using simmulated annealing heuristic.");
  try
  {
    options
      .add_options()
      ("h, help", "Show help")
      ("f, file", "Input file to read instance (ignore to read from stdin)", cxxopts::value<string>(), "file")
      ("s, start", "Start seed to use. Default is 0", cxxopts::value<int>(), "N")
      ("e, end", "End seed to use. Default is 10", cxxopts::value<int>(), "M")
      ("use-hybrid", "Use hybrid sweep calculation")
      ("skip-sweep", "Skip final sweep calculation");

    auto result = options.parse(argc, argv);
    if (result.count("h"))
    {
      cout << options.help() << endl;
      return 0;
    }
    if (result.count("f"))
    {
      inFile.open(result["f"].as<string>());
      is = &inFile;
    }
    if (result.count("s"))
    {
      sseed = result["s"].as<int>();
    }
    if (result.count("e"))
    {
      eseed = result["e"].as<int>();
    }
    if (result.count("use-hybrid"))
    {
      hybrid = true;
    }
    if (result.count("final-sweep"))
    {
      sweep = false;
    }
  } catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return 1;
  }

  // NUMBER OF CITIES
  const int n = N;
  Annealing annealing(n);

  // --DATABASE READ SECTION--
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;
  rc = sqlite3_open(DBPATH, &db);
  char sql[200];

  if (rc)
  {
    cout << " Can't open database: " << sqlite3_errmsg(db) << endl;
    return 0;
  }

  sprintf(sql, "SELECT * \
        FROM (SELECT * FROM cities LIMIT %d) A \
        LEFT JOIN (SELECT * FROM connections) B \
        ON A.id=B.id_city_1 AND B.id_city_2<=%d;",
          n, n);

  rc = sqlite3_exec(db, sql, callback, &annealing, &zErrMsg);

  if (rc != SQLITE_OK)
  {
    cout << "SQL error: " << zErrMsg << endl;
    sqlite3_free(zErrMsg);
  }
  sqlite3_close(db);
  // --END DATABASE READ SECTION--

  string s;
  vector<int> S;

  // READ TEST CASES FROM INPUT FILE (COMMA SEPARATED).
  while (getline(*is, s))
  {
    // --PARSE INPUT--
    stringstream ss;
    ss << s;
    while (getline(ss, s, ','))
    {
      S.push_back(stoi(s) - 1);
    }
    // --END PARSE INPUT--
    // --SIMULATED ANNEALING--
    for (int seed = sseed; seed < eseed; seed++)
    {
      auto start = chrono::steady_clock::now();
      // params: (seed, uniform_int_generator_max)
      annealing.setRandomEngine(seed, S.size());
      // params: (initial_instance, hybrid_sweep?, final_sweep?)
      pair<vector<int>, double> res = annealing.computeSolution(S, hybrid, sweep);

      printf("\nSeed: %d\n", seed);
      for (auto i : res.first)
      {
        printf("%d, ", i+1);
      }
      printf("\nEvaluation: %2.9f\n", res.second);
      auto end = chrono::steady_clock::now();
      printf("Elapsed time: %lld\n", chrono::duration_cast<chrono::seconds>(end - start).count());
    }
    S.clear();
    // --END SIMULATED ANNEALING--
  }
  auto gend = chrono::steady_clock::now();
  printf("Total elapsed time: %lld\n", chrono::duration_cast<chrono::seconds>(gend - gstart).count());
  return 0;
}

static int callback(void *annealing, int argc, char **argv, char **azColName)
{
  Annealing *a = static_cast<Annealing *>(annealing);
  if (argv[7] != NULL)
    a->validEdge(stoi(argv[0]) - 1, stoi(argv[7]) - 1, stod(argv[8]));
  a->addCity(stoi(argv[0]) - 1, {stod(argv[4]), stod(argv[5])});
  return 0;
}
