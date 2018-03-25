#pragma once
#include <iostream>
#include <ilcplex/ilocplex.h>

using namespace std;

class MultiwayCut {
private:
	/* the # of vertices of input graph G */
	int n_vertices;

	/* the # of terminals, i.e. dimension of our problem space */
	int n_terminals;

	/* weight of edges: w[n_vertices][n_vertices] */
	double **weight_matrix;

	/* output of LP-solver, i.e. vertices on the simplex: u[n_vertices][n_terminals] */
	double **simplex_vertices;

	/* terminal array: T[n_terminals] */
	int *terminals;

	/* input graph topology: G[n_vertices][n_vertices] */
	bool **edge_matrix;

public:
	/* Constructor */
	MultiwayCut(void);

	/* LP solver function */
	double LP_solver(void);

	/* rounding algorithm [exponential-clock ver] */
	double rounding_alg_exp(void);

	/* rounding algorithm [simplex distortion ver] */
	double rounding_alg_dist(void);

	/* wrapper rounding algorithm */
	double rounding_alg(void);
};