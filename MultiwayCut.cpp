#include "MultiwayCut.h"

MultiwayCut::MultiwayCut(void)
{
	srand(unsigned(time(0)));
	n_vertices = (rand() % MAX_N_VERTICES) + 1; // random 1<= x <= 1000
	n_terminals = (rand() % n_vertices) + 1; // random 1<= x <= vertices
	cout << "n_vertices = " << n_vertices << endl;
	cout << "n_terminals = " << n_terminals << endl;
	weight_matrix = new double*[n_vertices]; // random,  upper triangle matrix and diagonal elements is zero
	for (int i = 0; i < n_vertices; i++) {
		weight_matrix[i] = new double[n_vertices];
		for (int j = 0; j <= i; j++) {
			weight_matrix[i][j] = 0;
		}
		for (int j = 1; j < n_vertices - i; j++) {
			weight_matrix[i][j + i] = ((double)rand() / RAND_MAX) * WEIGHT_MAX;
		}
	}
	cout << "---- weight matrix ----" << endl;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = 0; j < n_vertices; j++) {
			cout << setw(7) << weight_matrix[i][j] << " ";
		}
		cout << endl;
	}
	simplex_vertices = new double*[n_vertices]; // initialize 0
	for (int i = 0; i < n_vertices; i++) {
		simplex_vertices[i] = new double[n_terminals];
		for (int j = 0; j < n_terminals; j++) {
			simplex_vertices[i][j] = 0;
		}
	}
	cout << "terminals index : ";
	terminals = new int[n_terminals]; // indices of vertices
	terminal_random_choice(); // assign indices of vertices to terminals
	for (int i = 0; i < n_terminals; i++) {
		cout << terminals[i] << " ";
	}
	cout << endl;
	edge_matrix = new bool*[n_vertices]; // random true or false
	//cout << "---- edge matrix----" << endl;
	for (int i = 0; i < n_vertices;){
		edge_matrix[i] = new bool[n_vertices];
		//cout << i << ": ";
		for (int j = 0; j <= i; j++) {
			edge_matrix[i][j] = false;
			//cout << setw(7) << edge_matrix[i][j] << " ";
		}
		for (int j = 1; j < n_vertices - i; j++) {
			if (rand() % 10 == 0) {
				edge_matrix[i][j+i] = true;
				//cout << setw(7) << edge_matrix[i][j+i] << " ";
			}
			else {
				edge_matrix[i][j+i] = false;
				//cout << setw(7) << edge_matrix[i][j + i] << " ";
			}
		}
		if (check_vertex_isolated(i) == false || n_vertices == 1 || i == n_vertices - 1) {
			i++;
		}
		else{
			delete edge_matrix[i];
		}
		//cout << endl;		
	}
	/* avoid the last vertex disconnected */ 
	while (check_vertex_isolated(n_vertices - 1) && n_vertices != 1) {
		for (int j = 0; j < n_vertices - 1; j++) {
			if (rand() % 10 == 0) {
				edge_matrix[j][n_vertices - 1] = true;
			}
			else {
				edge_matrix[j][n_vertices - 1] = false;
			}
		}
	}
	cout << "---- edge matrix----" << endl;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = 0; j < n_vertices; j++) {
			cout << setw(7) << edge_matrix[i][j] << " ";
		}
		cout <<" ( "<< check_vertex_isolated(i) << " <==  1-isolated, 0-connected)" << endl;
	}

	for (int i = 0; i < n_vertices; i++) {
		for (int j = 0; j < n_vertices; j++) {
			weight_matrix[i][j] *= edge_matrix[i][j];
		}
	}
	cout << "---- modified weight matrix ---" << endl;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = 0; j < n_vertices; j++) {
			cout << setw(7) << weight_matrix[i][j] << " ";
		}
		cout << endl;
	}
}

double MultiwayCut::LP_solver(void)
{
	return 0;
}

double MultiwayCut::rounding_alg_exp(void)
{
	return 0;
}

double MultiwayCut::rounding_alg_dist(void)
{
	return 0;
}

double MultiwayCut::rounding_alg(void)
{
	return 0;
}

void MultiwayCut::terminal_random_choice() {
	int count = 0;
	bool *check = new bool[n_vertices];
	for (int i = 0; i < n_vertices; i++) {
		check[i] = false;
	}
	for (int i = 0; i < n_terminals; ) {
		int n = rand() % n_vertices;
		if (check[n] == false) {
			check[n] = true;
			i++;
		}
	}
	for (int i = 0, j = 0; i < n_vertices; i++) {
		if (check[i] == true) {
			terminals[j] = i;
			j++;
		}
	}
}

bool MultiwayCut::check_vertex_isolated(int k) {
	int i = k + 1;
	for (; i < n_vertices; i++) {
		if (edge_matrix[k][i] == 1)
			break;
	}
	if (i == n_vertices) {
		int j = 0;
		for (; j < k; j++) {
			if (edge_matrix[j][k] == 1)
				break;
		}
		if (j == k) {
			return true;
		}
	}
	return false;
}