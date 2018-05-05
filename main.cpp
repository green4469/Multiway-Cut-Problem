#include "MultiwayCut.h"

int main(int argc, char* argv[]) {
	/* LP saves the relaxed solution, RS saves the rounded solution */
	double LP;
	double RS;

	srand((unsigned)time(NULL));

	/* If no argument, return with error message */
	if (argc == 1) {
		cout << "You need to give a parameter!" << endl;
		cout << "For example, './MultiwayCut 0001.txt'" << endl;
		return -1;
	}

	if (argc == 2) {
		/* output files */
		ofstream fout_simplex("MCP_OUT\\MCP_OUT_SIMPLEX.TXT");
		ofstream fout_summary("MCP_OUT\\MCP_OUT.TXT", ofstream::out | ofstream::app);
		ofstream fout_edge_cut("MCP_OUT\\MCP_OUT_EDGE_CUT.TXT", ofstream::out | ofstream::app);

		MultiwayCut *MC = new MultiwayCut(argc, argv);
		LP = MC->LP_solver();
		RS = MC->rounding_alg();
		cout << "relaxed solution: " << LP << endl;
		cout << "rounded solution: " << RS << endl;

		cout << argv[1] << endl;
		for (int i = 0; i < MC->n_vertices; i++) {
			cout << i << " = (";
			for (int j = 0; j < MC->n_terminals - 1; j++) {
				cout << MC->simplex_vertices[i][j] << ", ";
			}
			cout << MC->simplex_vertices[i][MC->n_terminals - 1] << ')' << endl;
		}
		cout << endl;

		/* simplex output */
		/*
		fout_simplex << argv[1] << endl;
		for (int i = 0; i < MC->n_vertices; i++) {
			fout_simplex << i << " = (";
			for (int j = 0; j < MC->n_terminals - 1; j++) {
				fout_simplex << MC->simplex_vertices[i][j] << ", ";
			}
			fout_simplex << MC->simplex_vertices[i][MC->n_terminals-1] << ')' << endl;
		}
		fout_simplex << endl;
		*/
		/* summary output */
		//fout_summary << file_num << "," << MC->n_vertices << "," << MC->n_terminals << "," << RS/LP << endl;

		/* edge_cut output */

		delete MC;
	}

	return 0;
}