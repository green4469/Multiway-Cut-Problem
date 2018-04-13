#include "MultiwayCut.h"

string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
);

int main(void) {
	double LP;
	double RS;
	int iteration = 0;

	ofstream fout("output.txt");
	int num_of_files;
	cout << "The number of files: ";
	cin >> num_of_files;
	cout << "The start number of files: ";
	int start_num;
	cin >> start_num;
	/* loop until different between relaxed solution & rounded solution	*/
	for(int i = 0; i < num_of_files; i++){
		string in_file = "MCP_IN\\MCP_IN_";
		string number;
		if (i / 10 == 0) {
			number = "000" + to_string(i);
		}
		else if (i / 100 == 0) {
			number = "00" + to_string(i);
		}
		else if (i / 1000 == 0) {
			number = "0" + to_string(i);
		}
		else {
			number = to_string(i);
		}
		in_file.append(number).append(".TXT");
		string out_file = "MCP_OUT\\MCP_OUT.TXT";
		//out_file = in_file;
		//out_file = replace_all(out_file, "IN", "OUT");
		ofstream fout(out_file,ofstream::out | ofstream::app);
		srand((unsigned)time(NULL) + (unsigned)iteration * 10);
		cout << ++iteration << "th case" << endl << endl;
		MultiwayCut *a;
		if(in_file.size() > 0)
			a = new MultiwayCut(in_file);
		else
			a = new MultiwayCut();
		LP = a->LP_solver();
		RS = a->rounding_alg();

		double duality_gap;
		if (CompareDoubleUlps(LP, 0) == 0 && CompareDoubleUlps(RS, 0) == 0)
			duality_gap = 0;
		else
			duality_gap = RS / LP;
		fout << a->n_vertices << "," << a->n_terminals << "," << duality_gap << endl;

		/*
		if (CompareDoubleUlps(LP, RS) != 0) {
			fout << "LP's objective value : " << LP << endl;
			fout << "Rounded objective value : " << RS << endl;
			fout << "Weight Matrix" << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					fout << a->weight_matrix[i][j] << '\t';
				}
				fout << endl;
			}

			fout << "Edge Matrix" << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					fout << a->edge_matrix[i][j] << '\t';
				}
				fout << endl;
			}

			fout << "Removed edge " << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					fout << a->removed_edge[i][j] << '\t';
				}
				fout << endl;
			}

			fout << "assigned terminal (lu) " << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				fout << a->assigned_terminal[i] << '\t';
			}
			fout << endl;

			fout << "Terminals" << endl;
			for (int i = 0; i < a->n_terminals; i++) {
				fout << a->terminals[i] << '\t';
			}
			fout << endl;
		}
		*/

		cout << "relaxed_solution: " << LP << endl;
		//cout << "optimal_solution: " << a.get_optimal_solution() << endl;
		cout << "rounded_solution: " << RS << endl;
		delete a;
	}

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

string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
) {

	std::string result = message;
	std::string::size_type pos = 0;
	std::string::size_type offset = 0;

	while ((pos = result.find(pattern, offset)) != std::string::npos)
	{
		result.replace(result.begin() + pos, result.begin() + pos + pattern.size(), replace);
		offset = pos + replace.size();
	}

	return result;
}

