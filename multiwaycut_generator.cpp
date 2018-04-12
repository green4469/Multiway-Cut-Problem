#include <iostream>
#include <fstream>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <iomanip>
#include <string>
#include <vector>
#include <windows.h>
#define MAX_N_VERTICES 10
#define WEIGHT_MAX 100
using namespace std;

void terminal_random_choice();
bool check_vertex_isolated(int k);

/* the # of vertices of input graph G */
int n_vertices = 0;
/* the # of terminals, i.e. dimension of our problem space */
int n_terminals;
/* weight of edges: w[n_vertices][n_vertices] */
double **weight_matrix;
/* terminal array: T[n_terminals] */
int *terminals;
/* input graph topology: G[n_vertices][n_vertices] */
bool **edge_matrix;
/* grouping array for making one graph */
int* group = NULL;
void grouping(vector<int>* group, int n_vertices);
void make_one_graph(vector<int>* group);
string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
);
bool check_array_all_true(bool arr[], int array_size);
bool check_num_in_group(vector<int> group, int n);
void MakeDirectory(string full_path);
int main(void) {
	srand(unsigned(NULL));
	int num_of_files;
	cout << "파일 몇개 만들까?: ";
	cin >> num_of_files;
	cout << "파일 시작 번호: ";
	int start_num;
	cin >> start_num;
	cout << "vertex 최소 갯수: ";
	int vertex_min;
	cin >> vertex_min;
	cout << "vertex 최대 갯수: ";
	int vertex_max;
	cin >> vertex_max;
	for (int i = start_num; i < start_num + num_of_files; i++) {
		srand(unsigned(NULL));
		string out_file = "MCP_IN\\MCP_IN_";
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
		out_file.append(number).append(".TXT");
		
		n_vertices = (rand() % (vertex_max - vertex_min + 1 )) + vertex_min; // random 1<= x <= 1000
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
		cout << "terminals index : ";

		terminals = new int[n_terminals]; // indices of vertices
		terminal_random_choice(); // assign indices of vertices to terminals
		for (int i = 0; i < n_terminals; i++) {
			cout << terminals[i] << " ";
		}
		cout << endl;

		cout << "---- edge matrix----" << endl;

		edge_matrix = new bool*[n_vertices]; // random true or false
		for (int i = 0; i < n_vertices; i++) {
			edge_matrix[i] = new bool[n_vertices];
			cout << i << ": ";
			for (int j = 0; j <= i; j++) {
				edge_matrix[i][j] = false;
				cout << setw(7) << edge_matrix[i][j] << " ";
			}
			for (int j = 1; j < n_vertices - i; j++) {
				if (rand() % 10 == 0) {
					edge_matrix[i][j + i] = true;
					cout << setw(7) << edge_matrix[i][j + i] << " ";
				}
				else {
					edge_matrix[i][j + i] = false;
					cout << setw(7) << edge_matrix[i][j + i] << " ";
				}
			}
			cout << endl;
		}

		cout << "---- edge matrix----" << endl;
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				cout << setw(7) << edge_matrix[i][j] << " ";
			}

			cout << " ( " << check_vertex_isolated(i) << " <==  1-isolated, 0-connected)" << endl;

		}
		vector<int>* group = new vector<int>[n_vertices];
		for (int i = 0; i < group[0].size(); i++) {
			cout << group[0][i] << " ";
		}
		cout << endl;
		grouping(group, n_vertices);
		for (int i = 0; i < group[0].size(); i++) {
			cout << group[0][i] << " ";
		}
		cout << endl;
		make_one_graph(group);
		for (int i = 0; i < group[0].size(); i++) {
			cout << group[0][i] << " ";
		}
		cout << endl;

		cout << "---- modified edge matrix----" << endl;
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				cout << setw(7) << edge_matrix[i][j] << " ";
			}
			cout << " ( " << check_vertex_isolated(i) << " <==  1-isolated, 0-connected)" << endl;
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

		MakeDirectory(out_file);
		ofstream out(out_file);
		out << n_vertices << endl;
		out << n_terminals << endl;
		for (int i = 0; i < n_terminals; i++) {
			out << terminals[i] << " ";
		}
		out << endl;
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				if (weight_matrix[i][j] != 0) {
					out << i << " " << j << " " << weight_matrix[i][j] << endl;
				}
			}
		}
		delete[] group;
		Sleep(100);
	}
	return 0;
}

void terminal_random_choice() {
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
	delete check;
}

bool check_vertex_isolated(int k) {
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


void grouping(vector<int>* group, int n_vertices) {
	bool *check;
	int group_number = 0;
	check = new bool[n_vertices];
	for (int i = 0; i < n_vertices; i++) {
		check[i] = false;
	}
	check[0] = true;
	vector<int> one_group_queue;
	for (int i = 0; i < n_vertices; i++) {
		if (check_array_all_true(check, n_vertices))
			break;
		group[group_number].push_back(i);
		for (int j = 0; j < n_vertices; j++) {
			if (check[j] == false && edge_matrix[i][j] == 1) {
				group[group_number].push_back(j);
				check[j] = true;
				one_group_queue.push_back(j);
			}
		}
		while (one_group_queue.size() != 0) {
			int vertex = one_group_queue.back();
			one_group_queue.pop_back();
			for (int j = 0; j < n_vertices; j++) {
				if (check[j] == false && edge_matrix[vertex][j] == 1) {
					group[group_number].push_back(j);
					check[j] = true;
					one_group_queue.push_back(j);
				}
			}
		}
		group_number++;
	}
	delete check;
}

bool check_array_all_true(bool arr[], int array_size) {
	for (int i = 0; i < array_size; i++) {
		if (arr[i] == false) {
			return false;
		}
	}
	return true;
}

void make_one_graph(vector<int>* group) {
	srand(time(NULL));
	for (int i = 1; i < n_vertices; i++) {
		bool check_connect = false;
		if (group[i].size() != 0) {
			do {
				for (int j = 0; j < group[0].size(); j++) {
					for (int k = 0; k < group[i].size(); k++) {
						if (rand() % 10 == 0) {
							if (group[0][j] < group[i][k]) {
								edge_matrix[group[0][j]][group[i][k]] = 1;
							}
							else {
								edge_matrix[group[i][k]][group[0][j]] = 1;
							}
							check_connect = true;
						}
					}
				}
			} while (check_connect == false);
			for (int j = 0; j < group[i].size(); j++) {
				if(!check_num_in_group(group[0], group[i][j]))
					group[0].push_back(group[i][j]);
			}
		}
	}
	return;
}

bool check_num_in_group(vector<int> group, int n){
	bool check = false;
	for (int i = 0; i < group.size(); i++) {
		if (group[i] == n) {
			return true;
		}
	}
	return false;
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

void MakeDirectory(string full_path)
{
	string strDirPath;
	string strTempDir;
	int indexOf = 0;
	while (true) {
		indexOf = strDirPath.find("\\");
		strTempDir += strDirPath.substr(0, indexOf) + "\\";
		CreateDirectory(strTempDir.c_str(),NULL);
		strDirPath = strDirPath.substr(indexOf + 1, strDirPath.length());
		if (indexOf < 0) break;
	}
}