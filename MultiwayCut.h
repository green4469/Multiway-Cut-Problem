#pragma once
#include <iostream>
#include <cstdlib>
#include <ilcplex/ilocplex.h>
#include <random>
#include <string>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>
#include <numeric>		// std::iota
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /* time */
#include <fstream>

#define MAX_N_VERTICES 100
#define WEIGHT_MAX 100
using namespace std;

class MultiwayCut {
public:
//private:
	/* the # of vertices of input graph G */
	int n_vertices;

	/* the # of terminals, i.e. dimension of our problem space */
	int n_terminals;

	/* weight of edges: w[n_vertices][n_vertices] */
	double **weight_matrix;

	/* output of LP-solver, i.e. vertices on the simplex: u[n_vertices][n_terminals] */
	double **simplex_vertices;

	/* optimal solution [n_vertices][n_terminals]*/
	double **optimal_solution;

	/* terminal array: T[n_terminals] */
	int *terminals;

	/* input graph topology: G[n_vertices][n_vertices] */
	bool **edge_matrix;

	/* denote that each vertex assigned which terminal. in pseudo code, l(u). l[n_vertices] */
	int *assigned_terminal;

	/* upper triangle matrix with diagonal elements '0' */
	bool **removed_edge;


//public:
	/* Constructor */
	MultiwayCut(int argc, char *argv[]);

	/* Deconstructor */
	~MultiwayCut(void);

	/* LP solver function */
	double LP_solver(void);

	/* rounding algorithm [exponential-clock ver] */
	double rounding_alg_exp(void);

	/* rounding algorithm [simplex distortion ver] */
	double rounding_alg_dist(void);

	/* wrapper rounding algorithm */
	double rounding_alg(void);

	/* choose terminals from vertices */
	void terminal_random_choice();

	/* if k'th vertex has no edges return true, else return false */
	bool check_vertex_isolated(int k);

	/* Post-processing function*/
	double post_process(void);
};

int CompareDoubleUlps(double x, double y, int ulpsTolerance = 1000000);
