#include "MultiwayCut.h"

int main(int argc, char* argv[]) {
	/* LP saves the relaxed solution, RS saves the rounded solution */
	double LP;
	double RS;

	if (argc == 2) {
		ofstream of("output.txt");

		MultiwayCut *MC = new MultiwayCut(argc, argv);

		LP = MC->LP_solver();
		RS = MC->rounding_alg();
		cout << "relaxed solution: " << LP << endl;
		cout << "rounded solution: " << RS << endl;

		of << "Input graph" << endl;
		for (int i = 0; i < MC->n_vertices; i++) {
			for (int j = 0; j < MC->n_vertices; j++) {
				of << MC->weight_matrix[i][j] << '\t';
			}
			of << endl;
		}

		of << endl;

		of << "Coordinate values of vertices on simplex" << endl;
		for (int i = 0; i < MC->n_vertices; i++) {
			of << i << " = (";
			for (int j = 0; j < MC->n_terminals - 1; j++) {
				of << MC->simplex_vertices[i][j] << ", ";
			}
			of << MC->simplex_vertices[i][MC->n_terminals-1] << ')' << endl;
		}

		of << endl;

		of << "Removed edges" << endl;
		for (int i = 0; i < MC->n_vertices; i++) {
			for (int j = 0; j < MC->n_vertices; j++) {
				of << MC->removed_edge[i][j] << ' ';
			}
			of << endl;
		}
	}

	/*
	MultiwayCut a;
	LP = a.LP_solver();
	RS = a.rounding_alg();
	cout << "relaxed solution: " << LP << endl;
	cout << "rounded solution: " << RS << endl;
	*/


	/* loop until different between relaxed solution & rounded solution	(random-generation)
	ofstream of("output.txt");
	int iteration = 0;
	do {
		srand((unsigned)time(NULL) + (unsigned)iteration * 10);
		cout << ++iteration << "th case" << endl << endl;
		MultiwayCut *a = new MultiwayCut;
		LP = a->LP_solver();
		RS = a->rounding_alg();

		if (CompareDoubleUlps(LP, RS) != 0) {
			of << "LP's objective value : " << LP << endl;
			of << "Rounded objective value : " << RS << endl;
			of << "Weight Matrix" << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					of << a->weight_matrix[i][j] << '\t';
				}
				of << endl;
			}

			of << "Edge Matrix" << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					of << a->edge_matrix[i][j] << '\t';
				}
				of << endl;
			}

			of << "Removed edge " << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					of << a->removed_edge[i][j] << '\t';
				}
				of << endl;
			}

			of << "assigned terminal (lu) " << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				of << a->assigned_terminal[i] << '\t';
			}
			of << endl;

			of << "Terminals" << endl;
			for (int i = 0; i < a->n_terminals; i++) {
				of << a->terminals[i] << '\t';
			}
			of << endl;
		}

		cout << "relaxed_solution: " << LP << endl;
		//cout << "optimal_solution: " << a.get_optimal_solution() << endl;
		cout << "rounded_solution: " << RS << endl;
		delete a;
	} while (CompareDoubleUlps(LP,RS) == 0);

	cout << "Different!" << endl;
	*/

	return 0;
}