#include "MultiwayCut.h"

int main(void) {
	double LP;
	double RS;
	int iteration = 0;

	ofstream of("output.txt");

	/* loop until different between relaxed solution & rounded solution	*/
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

	/*
	MultiwayCut a;
	LP = a.LP_solver();
	RS = a.rounding_alg();
	cout << "relaxed solution: " << LP << endl;
	cout << "rounded solution: " << RS << endl;
	*/

	//cout << "Different!" << endl;
	return 0;
}