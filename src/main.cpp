#include <iostream>
#include <cmath>
#include <sstream>
#include <sqlite3.h>
#include "Graph.h"
#include "Annealing.h"

using namespace std;

// DB FILE PATH
#define DBPATH "db/tsp.db"

// NUMBER OF CITIES TO READ FROM DB
#define N 1092;

template <int n>
static int callback(void *z, int argc, char **argv, char **azColName)
{
  if (argv[7] != NULL)
    Recocido<n>::aristaValida(stoi(argv[0]) - 1, stoi(argv[7]) - 1, stod(argv[8]));
  Recocido<n>::addCity(stoi(argv[0]) - 1, {stod(argv[4]), stod(argv[5])});
  return 0;
}

int main()
{
  // NUMBER OF CITIES
  const int n = N;

  // --DATABASE READ SECTION--
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;
  rc = sqlite3_open(DBPATH, &db);
  char const *sql;

  if (rc)
  {
    cout << " Can't open database: " << sqlite3_errmsg(db) << endl;
    return 0;
  }
  else
  {
    cout << "Opened database successfully" << endl;
  }

  sql = "SELECT * \
        FROM (SELECT * FROM cities LIMIT 1092) A \
        LEFT JOIN (SELECT * FROM connections) B \
        ON A.id=B.id_city_1 AND B.id_city_2<=1092;";

  rc = sqlite3_exec(db, sql, callback<n>, 0, &zErrMsg);

  if (rc != SQLITE_OK)
  {
    cout << "SQL error: " << zErrMsg << endl;
    sqlite3_free(zErrMsg);
  }
  else
  {
    cout << "Operation done successfully" << endl;
  }
  sqlite3_close(db);

  string s;
  vector<int> S;
  int testcase = 0;

  // READ TEST CASES FROM STDIN (COMMA SEPARATED).
  while (getline(cin, s))
  {
    // --PARSE INPUT--
    stringstream ss;
    ss << s;
    while(getline(ss, s, ',')) {
      S.push_back(stoi(s)-1);
    }
    // --END PARSE INPUT--

    // --SIMULATED ANNEALING--
    double res = Recocido<n>::costFunction(S);
    printf("Evaluation#%d: %2.9f\n\n", ++testcase, res);
    S.clear();
    // --END SIMULATED ANNEALING--
  }

  return 0;
}
