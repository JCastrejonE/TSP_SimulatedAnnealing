## About
This project forms part of Heuristic combinatorial optimization course that forms part of Computer Science degree at **Facultad de Ciencias, Universidad Nacional Autónoma de México**. It tries to find a solution to Traveling salesman problem using simmulated annealing algorithm with some few variants.

For further understanding you can refer to the [PDF document in spanish](tex/Recocido_Simulado.pdf).

## Building and running

[Meson build system](https://mesonbuild.com/index.html) is used to compile this project. Running:
```
$ meson build && cd build
$ ninja
```
in the root directory should setup *build* directory with the executable inside.

Use `$ build/TSP_SimmulatedAnnealing -h` to know more about the available options.
```
$ ./TSP_SimmulatedAnnealing -h
TSP problem solution using simmulated annealing heuristic.
Usage:
  ./TSP_SimmulatedAnnealing [OPTION...]

  -h, --help              Show help
  -f, --file input_file   Input file to read instance (ignore to read from
                          stdin)
  -s, --start N           Start seed to use. Default is 0
  -e, --end M             End seed to use. Default is 10
  -p, --plot output_file  Create gnuplot script
      --use-hybrid        Use hybrid sweep calculation
      --skip-sweep        Skip final sweep calculation
```

## Plot option
This option generates a [gnuplot](http://gnuplot.sourceforge.net) script that shows an historic of each accepted solution in the algorithm. Use `$ gnuplot plot_file_name` to get the plotted graph.

## Available tune options
The simmulated annealing heuristic can be tuned to experiment different results.
The header file [include/Annealing.hpp](include/Annealing.hpp) contains the constants that can be modified.

Credits to [jarro2783](https://github.com/jarro2783/) for his [cxxopts](https://github.com/jarro2783/cxxopts) library to parse command line options.