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

	int file_num = -1;
	while(true){
		ifstream fin;
		string in_file = "MCP_IN\\MCP_IN_";
		file_num++;
		string number;
		if (file_num / 10 == 0) {
			number = "000" + to_string(file_num);
		}
		else if (file_num / 100 == 0) {
			number = "00" + to_string(file_num);
		}
		else if (file_num / 1000 == 0) {
			number = "0" + to_string(file_num);
		}
		else {
			number = to_string(file_num);
		}
		in_file.append(number).append(".TXT");
		fin.open(in_file, ifstream::in);
		if (fin.fail() == true)
			break;
	}

	/* loop until different between relaxed solution & rounded solution	*/
	do {
		string in_file = "MCP_IN\\MCP_IN_";
		string number;
		if (file_num / 10 == 0) {
			number = "000" + to_string(file_num);
		}
		else if (file_num / 100 == 0) {
			number = "00" + to_string(file_num);
		}
		else if (file_num / 1000 == 0) {
			number = "0" + to_string(file_num);
		}
		else {
			number = to_string(file_num);
		}
		in_file.append(number).append(".TXT");
		string out_file_summary = "MCP_OUT\\MCP_OUT.TXT";
		//string out_file = in_file;
		//out_file = replace_all(out_file, "IN", "OUT");
		string out_file = "MCP_OUT\\MCP_OUT_EDGE_CUT.TXT";
		srand((unsigned)time(NULL) + (unsigned)iteration * 10);
		cout << ++iteration << "th case" << endl << endl;
		MultiwayCut *a;
		a = new MultiwayCut();
		cout << ++iteration << "constructor" << endl << endl;
		LP = a->LP_solver();
		RS = a->rounding_alg();
		
		if (CompareDoubleUlps(LP, RS) != 0) {
			ofstream fout_result(out_file, ofstream::out | ofstream::app);
			ofstream fout_one(in_file);
			ofstream fout_summary(out_file_summary, ofstream::out | ofstream::app);
			double duality_gap;
			duality_gap = RS / LP;
			/* summary file */
			fout_summary << file_num << "," << a->n_vertices << "," << a->n_terminals << "," << duality_gap << endl;
			/* each input file */
			fout_one << a->n_vertices << endl;
			fout_one << a->n_terminals << endl;
			for (int i = 0; i < a->n_terminals; i++) {
				fout_one << a->terminals[i] << " ";
			}
			fout_one << endl;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					if (a->weight_matrix[i][j] != 0) {
						fout_one << i << " " << j << " " << a->weight_matrix[i][j] << endl;
					}
				}
			}
			/* result file */
			fout_result << file_num << "," << LP << "," << RS << "," << RS / LP;
			for (int i = 0; i < a->n_vertices; i++) {
				for (int j = 0; j < a->n_vertices; j++) {
					if (a->removed_edge[i][j] == true) {
						fout_result << "," << i << "-" << j;
					}
				}
			}
			fout_result << endl;
			file_num++;
		}

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
	} while (1);// CompareDoubleUlps(LP, RS) == 0);

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

