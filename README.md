# MAXE - a MOLP helper for Maximum Entropy Method

This version of [outer](https://github.com/lcsirmaz/outer) is a helper
program for determining consequences of the Maximum Entropy Method.
The main differences are:

* objectives must be non-negative
* all facets, vertices, and extremal directions of the objective polyhedron are determined
* an initial internal point must be specified in the VLP file

#### USAGE

The program is invoked as

    maxe [options] <vlp-file>

The only obligatory argument is the file name which contains the description
of the problem in vlp format. Accepted options and configuration parameters
are the same as for the `outer` program. Use

    maxe --help=vlp

for a description how the initial internal point is defined.

#### COMPILATION

The program uses glpk, the GNU Linear Program Kit, for solving scalar LP problems.
Changing to the directory `MAXE/src`, the following command compiles **maxe** 
linking the glpk routines:

    gcc -O3 -W -o maxe *.c -lm -lglpk

#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>


