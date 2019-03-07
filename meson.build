project(
  'TCP_SimmulatedAnnealing',
  'cpp', default_options : ['cpp_std=c++11'])

incdir = include_directories('include')

sqlite3 = dependency('sqlite3')

executable(
  'TCP_SimmulatedAnnealing',
  'src/main.cpp',
  include_directories : incdir,
  dependencies : sqlite3)