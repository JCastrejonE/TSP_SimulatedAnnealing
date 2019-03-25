project(
  'TSP_SimmulatedAnnealing',
  'cpp', default_options : ['cpp_std=c++11'])

incdir = include_directories('include')

sqlite3 = dependency('sqlite3')

src = ['src/Annealing.cpp', 'src/Graph.cpp', 'src/Main.cpp']

executable(
  'TSP_SimmulatedAnnealing',
  src,
  include_directories : incdir,
  dependencies : sqlite3)