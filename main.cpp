#include "MultiwayCut.h"

int main(int argc, char* argv[]) {
	/* LP saves the relaxed solution, RS saves the rounded solution */
	double LP;
	double RS;
	clock_t RT_LP, RT_RD;

	srand((unsigned)time(NULL));
	rand();

	/* If no argument, return with error message */
	if (argc == 1) {
		cout << "You need to give a parameter!" << endl;
		cout << "For example, './MultiwayCut 0001.txt'" << endl;
		return -1;
	}

	if (argc == 2) {
		/* output files */
		ofstream fout_simplex("MCP_OUT_SIMPLEX.TXT", ofstream::out | ofstream::app);
		ofstream fout_summary("MCP_OUT.TXT", ofstream::out | ofstream::app);
		//ofstream fout_edge_cut("MCP_OUT\\MCP_OUT_EDGE_CUT.TXT", ofstream::out | ofstream::app);

		string fnum;
		string fname = argv[1];
		fnum.append(fname, 9, 7);
		cout << "fnum: " << fnum << endl;

		MultiwayCut *MC = new MultiwayCut(argc, argv);
		RT_LP = clock();
		LP = MC->LP_solver();
		RT_LP = clock() - RT_LP;

		RT_RD = clock();
		RS = MC->rounding_alg();
		RT_RD = clock() - RT_RD;

		float RT_LP_S = ((float)RT_LP) / CLOCKS_PER_SEC;
		float RT_RD_S = ((float)RT_RD) / CLOCKS_PER_SEC;

		cout << "relaxed solution: " << LP << endl;
		cout << "rounded solution: " << RS << endl;
		cout << "IG: " << RS / LP << endl;

		fout_summary << fnum << ' ' << LP << ' ' << RS << ' ' << RS / LP << ' ' << RT_LP_S << ' ' << RT_RD_S << ' ' << RT_LP_S + RT_RD_S << endl;

		/* simplex output */
		//fout_simplex << argv[1] << endl;
		fout_simplex << LP << ' ' << RS  << ' ' << RS / LP << endl;

		for (int i = 0; i < MC->n_vertices; i++) {
			//fout_simplex << i << " = (";
			for (int j = 0; j < MC->n_terminals - 1; j++) {
				fout_simplex << MC->simplex_vertices[i][j] << ' ';
			}
			fout_simplex << MC->simplex_vertices[i][MC->n_terminals-1] << endl;
		}
		fout_simplex << endl;

		//delete MC;
	}

	return 0;
}