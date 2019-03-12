#include <iostream>
#include <cmath>
#include <sstream>
#include <chrono>
#include <sqlite3.h>
#include "Graph.hpp"
#include "Annealing.hpp"

using namespace std;

// DB FILE PATH
#define DBPATH "../db/tsp.db"

// NUMBER OF CITIES TO READ FROM DB
#define N 1092

// GETS EXECUTED FOR EACH ROW OF THE DB QUERY RESULT.
static int callback(void *, int, char **, char **);

int main()
{
  auto gstart = chrono::steady_clock::now();
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
  // int testcase = 0;

  // READ TEST CASES FROM STDIN (COMMA SEPARATED).
  while (getline(cin, s))
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
    for (int seed = 500; seed < 600; seed++)
    {
      auto start = chrono::steady_clock::now();
      annealing.setRandomEngine(seed, S.size());
      pair<vector<int>, double> res = annealing.computeSolution(S, true);

      printf("\nSeed: %d\n", seed);
      for (auto i : res.first)
      {
        printf("%d, ", i);
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
