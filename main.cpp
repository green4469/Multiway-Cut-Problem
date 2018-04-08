#include "MultiwayCut.h"
string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
);
int main(int argc, char* argv[]) { // argv : file name (ex: MCP_IN_0001.txt)
	double LP;
	double RS;
	string out_file;
	out_file = argv[1];
	out_file = replace_all(out_file,"IN", "OUT");
	ofstream out(out_file);
	/* loop until different between relaxed solution & rounded solution	*/
	//do {
		MultiwayCut a(argc, argv);
		LP = a.LP_solver();
		RS = a.rounding_alg();
		cout << "relaxed_solution: " << LP << endl;
		//cout << "optimal_solution: " << a.get_optimal_solution() << endl;
		cout << "rounded_solution: " << RS << endl;
	//} while (CompareDoubleUlps(LP,RS) == 0);
		out << LP << "," << RS;
		int n_vertices = a.get_n_vertices();
		bool** removed_edge = a.get_removed_edge();
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				if (removed_edge[i][j] == true) {
					out << "," << i << "-" << j;
				}
			}
		}
		out << endl;
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