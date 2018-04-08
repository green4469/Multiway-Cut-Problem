#include "MultiwayCut.h"

int CompareDoubleUlps(double x, double y, int ulpsTolerance)
{
	double diff = x - y;

	__int64 nx = *((__int64*)&x);
	__int64 ny = *((__int64*)&y);

	if ((nx & 0x8000000000000000) != (ny & 0x8000000000000000))
	{
		if (x == y)
			return 0;

		return (diff > 0) ? 1 : -1;
	}

	__int64 ulpsDiff = nx - ny;
	if ((ulpsDiff >= 0 ? ulpsDiff : -ulpsDiff) <= ulpsTolerance)
		return 0;

	return (diff > 0) ? 1 : -1;
}

MultiwayCut::MultiwayCut(int argc, char* argv[])
{
	optimal_solution = (double **)malloc(sizeof(double *) * n_vertices);
	for (int i = 0; i < n_vertices; i++) {
		optimal_solution[i] = (double *)malloc(sizeof(double) * n_terminals);
	}
	
	bool check = false;
	switch (argc) {
	case 1: // random sample
		srand(unsigned(time(0)));
		n_vertices = 6; // (rand() % MAX_N_VERTICES) + 1; // random 1<= x <= 1000
		n_terminals = 3; // (rand() % n_vertices) + 1; // random 1<= x <= vertices

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

		edge_matrix = new bool*[n_vertices]; // random true or false
		//cout << "---- edge matrix----" << endl;
		for (int i = 0; i < n_vertices;) {
			if (check == false)
				edge_matrix[i] = new bool[n_vertices];
			else
				check = false;
			///cout << i << ": ";
			for (int j = 0; j <= i; j++) {
				edge_matrix[i][j] = false;
				///cout << setw(7) << edge_matrix[i][j] << " ";
			}
			for (int j = 1; j < n_vertices - i; j++) {
				if (rand() % 2 == 0) {
					edge_matrix[i][j + i] = true;
					///cout << setw(7) << edge_matrix[i][j+i] << " ";
				}
				else {
					edge_matrix[i][j + i] = false;
					///cout << setw(7) << edge_matrix[i][j + i] << " ";
				}
			}
			if (check_vertex_isolated(i) == false || n_vertices == 1 || i == n_vertices - 1) {
				i++;
			}
			else {
				//delete edge_matrix[i];
				check = true;
			}
			///cout << endl;		
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

		simplex_vertices = new double*[n_vertices]; // initialize 0
		for (int i = 0; i < n_vertices; i++) {
			simplex_vertices[i] = new double[n_terminals];
			for (int j = 0; j < n_terminals; j++) {
				simplex_vertices[i][j] = 0;
			}
		}

		///cout << "---- edge matrix----" << endl;
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
		break;

	case 2: // file input
		ifstream in(argv[1]);
		//FILE * fp = fopen(argv[0], "r");
		in >> n_vertices;
		in >> n_terminals;
		cout << "n_vertices = " << n_vertices << endl;
		cout << "n_terminals = " << n_terminals << endl;
		int t;
		vector<int> temp;
		for (int i = 0; i < n_terminals; i++){
			in >> t;
			temp.push_back(t);
			//cout << "t: " << t << ", " << temp[i] << endl;
		}

		
		cout << "terminals index : ";
		terminals = new int[n_terminals]; // indices of vertices
		terminal_random_choice(); // assign indices of vertices to terminals
		for (int i = 0; i < n_terminals; i++) {
			terminals[i] = temp[i];
			cout << terminals[i] << " ";
		}
		cout << endl;
		//make weight_matrix filled 0
		weight_matrix = new double*[n_vertices]; // initialize 0
		for (int i = 0; i < n_vertices; i++) {
			weight_matrix[i] = new double[n_vertices];
			for (int j = 0; j < n_vertices; j++) {
				weight_matrix[i][j] = 0;
			}
		}
		int v1, v2;
		double w;
		
		while (in){
			in >> v1 >> v2 >> w;
			if (v1 < v2) {
				weight_matrix[v1][v2] = w;
			}
			else {
				weight_matrix[v2][v1] = w;
			}
		}
		///cout << "---- weight matrix ----" << endl;
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				cout << setw(7) << weight_matrix[i][j] << " ";
			}
			cout << endl;
		}

		terminals = new int[n_terminals];
		for (int i = 0; i < n_terminals; i++) {
			terminals[i] = temp[i];
		}

		simplex_vertices = new double*[n_vertices]; // initialize 0
		for (int i = 0; i < n_vertices; i++) {
			simplex_vertices[i] = new double[n_terminals];
			for (int j = 0; j < n_terminals; j++) {
				simplex_vertices[i][j] = 0;
			}
		}

		edge_matrix = new bool*[n_vertices]; // random true or false
		cout << "---- edge matrix----" << endl;
		for (int i = 0; i < n_vertices; i++) {
			edge_matrix[i] = new bool[n_vertices];
			//cout << i << ": ";
			for (int j = 0; j < n_vertices; j++) {
				if (weight_matrix[i][j] != 0)
					edge_matrix[i][j] = true;
				else
					edge_matrix[i][j] = false;
				//cout << setw(7) << edge_matrix[i][j] << " ";
			}
			cout << endl;
		}
		for (int i = 0; i < n_vertices; i++) {
			for (int j = 0; j < n_vertices; j++) {
				cout << setw(7) << edge_matrix[i][j] << " ";
			}
			cout << " ( " << check_vertex_isolated(i) << " <==  1-isolated, 0-connected)" << endl;
		}
	}
}

double MultiwayCut::get_optimal_solution(void) {
	IloEnv env;

	/* initialize mappings */
	IloNumVar **l = new IloNumVar*[n_vertices];
	for (int u = 0; u < n_vertices; ++u) {
		l[u] = new IloNumVar[n_terminals];
		for (int t = 0; t < n_terminals; ++t) {
			l[u][t] = IloNumVar(env, 0, 1, IloNumVar::Int);

		}
	}

	/* initialize L1 norms */
	IloNumVar ***z = new IloNumVar**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		z[i] = new IloNumVar*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			z[i][j] = new IloNumVar[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				z[i][j][k] = IloNumVar(env, -IloInfinity, IloInfinity);
			}
		}
	}


	/* set ranges of vertices */
	IloRange *vertexRanges = new IloRange[n_vertices];
	IloExpr *sumVertexComponents = new IloExpr[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		sumVertexComponents[i] = IloExpr(env);
		for (int j = 0; j < n_terminals; ++j)
			sumVertexComponents[i] += l[i][j];
		vertexRanges[i] = (sumVertexComponents[i] == 1.0);
	}

	/* set ranges of L1 norms */
	IloRange ***l1Ranges1 = new IloRange**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		l1Ranges1[i] = new IloRange*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			l1Ranges1[i][j] = new IloRange[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				l1Ranges1[i][j][k] = IloRange(env, 0, IloInfinity);
				l1Ranges1[i][j][k].setLinearCoef(z[i][j][k], 1);
				l1Ranges1[i][j][k].setLinearCoef(l[i][k], -1);
				l1Ranges1[i][j][k].setLinearCoef(l[j][k], 1);
			}
		}
	}
	IloRange ***l1Ranges2 = new IloRange**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		l1Ranges2[i] = new IloRange*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			l1Ranges2[i][j] = new IloRange[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				l1Ranges2[i][j][k] = IloRange(env, 0, IloInfinity);
				l1Ranges2[i][j][k].setLinearCoef(z[i][j][k], 1);
				l1Ranges2[i][j][k].setLinearCoef(l[i][k], 1);
				l1Ranges2[i][j][k].setLinearCoef(l[j][k], -1);
			}
		}
	}

	IloRange **connectionConstraint = new IloRange*[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		connectionConstraint[i] = new IloRange[n_terminals];
		for (int j = 0; j < n_terminals; ++j) {
			if (this->edge_matrix[i][terminals[j]] == false && this->edge_matrix[terminals[j]][i] == false) {
				connectionConstraint[i][j] = (l[i][j] == 0.0);
			}
			else {
				connectionConstraint[i][j] = IloRange(env, -IloInfinity, IloInfinity);
				connectionConstraint[i][j].setLinearCoef(l[i][j], 1);
			}
		}
	}


	/* objective function */
	IloObjective obj = IloMinimize(env, 0);
	for (int i = 0; i < n_vertices; ++i) {
		for (int j = 0; j < n_vertices; ++j) {
			for (int k = 0; k < n_terminals; ++k) {
				if (this->edge_matrix[i][j] == true)
					obj.setLinearCoef(z[i][j][k], 0.5 * this->weight_matrix[i][j]);
			}
		}
	}
	/* compile the model */
	IloModel model(env);
	for (int i = 0; i < n_vertices; ++i) {
		model.add(vertexRanges[i]);
		for (int j = 0; j < n_vertices; ++j) {
			for (int k = 0; k < n_terminals; ++k) {
				if (this->edge_matrix[i][j] == true) {
					model.add(l1Ranges1[i][j][k]);
					model.add(l1Ranges2[i][j][k]);
				}
			}
		}
	}
	for (int i = 0; i < n_vertices; ++i) {
		for (int j = 0; j < n_terminals; ++j) {
			if (this->edge_matrix[i][terminals[j]] == 0 && this->edge_matrix[terminals[j]][i] == 0) {
				model.add(connectionConstraint[i][j]);
			}
		}
	}
	model.add(obj);

	/* solve the model */
	IloCplex solver(model);
	try {
		solver.solve();
	}
	catch (IloException &ex) {
		cerr << ex << endl;
		return -1;
	}
	int stat = solver.getStatus();
	if (stat != 1) {
		cout << "solver status : " << stat << endl;
		return -1;
	}
	cout << solver.getObjValue() << endl;
	/* save results*/
	for (int i = 0; i < n_vertices; ++i) {
		cout << "vertex " << i << " : ";
		for (int j = 0; j < n_terminals; ++j) {
			this->optimal_solution[i][j] = solver.getValue(l[i][j]);
			cout << this->optimal_solution[i][j] << ' ';
		}
		cout << endl;
	}
	if (CompareDoubleUlps(solver.getObjValue(), 0.0) <=0) {
		return 0.0;
	}
	else
		return solver.getObjValue();

}

double MultiwayCut::LP_solver(void)
{
	IloEnv env;
	/* initialize vertices */
	IloNumVar **u = new IloNumVar*[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		u[i] = new IloNumVar[n_terminals];
		for (int j = 0; j < n_terminals; ++j) {
			u[i][j] = IloNumVar(env, 0, IloInfinity);
		}
	}
	/* initialize L1 norms */
	IloNumVar ***z = new IloNumVar**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		z[i] = new IloNumVar*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			z[i][j] = new IloNumVar[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				z[i][j][k] = IloNumVar(env, -IloInfinity, IloInfinity);
			}
		}
	}
	/* set ranges of terminals */
	IloRange **terminalRanges = new IloRange*[n_terminals];
	for (int i = 0; i < n_terminals; ++i) {
		terminalRanges[i] = new IloRange[n_terminals];
		for (int j = 0; j < n_terminals; ++j) {
			if (i == j)
				terminalRanges[i][j] = (u[terminals[i]][j] == 1.0);
			else
				terminalRanges[i][j] = (u[terminals[i]][j] == 0.0);
		}
	}
	/* set ranges of vertices */
	IloRange *vertexRanges = new IloRange[n_vertices];
	IloExpr *sumVertexComponents = new IloExpr[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		sumVertexComponents[i] = IloExpr(env);
		for (int j = 0; j < n_terminals; ++j)
			sumVertexComponents[i] += u[i][j];
		vertexRanges[i] = (sumVertexComponents[i] == 1.0);
	}
	/* set ranges of L1 norms */
	IloRange ***l1Ranges1 = new IloRange**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		l1Ranges1[i] = new IloRange*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			l1Ranges1[i][j] = new IloRange[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				l1Ranges1[i][j][k] = IloRange(env, 0, IloInfinity);
				l1Ranges1[i][j][k].setLinearCoef(z[i][j][k], 1);
				l1Ranges1[i][j][k].setLinearCoef(u[i][k], -1);
				l1Ranges1[i][j][k].setLinearCoef(u[j][k], 1);
			}
		}
	}
	IloRange ***l1Ranges2 = new IloRange**[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		l1Ranges2[i] = new IloRange*[n_vertices];
		for (int j = 0; j < n_vertices; ++j) {
			l1Ranges2[i][j] = new IloRange[n_terminals];
			for (int k = 0; k < n_terminals; ++k) {
				l1Ranges2[i][j][k] = IloRange(env, 0, IloInfinity);
				l1Ranges2[i][j][k].setLinearCoef(z[i][j][k], 1);
				l1Ranges2[i][j][k].setLinearCoef(u[i][k], 1);
				l1Ranges2[i][j][k].setLinearCoef(u[j][k], -1);
			}
		}
	}


	/* objective function */
	IloObjective obj = IloMinimize(env, 0);
	for (int i = 0; i < n_vertices; ++i) {
		for (int j = 0; j < n_vertices; ++j) {
			for (int k = 0; k < n_terminals; ++k) {
				if (this->edge_matrix[i][j] == true)
					obj.setLinearCoef(z[i][j][k], 0.5 * this->weight_matrix[i][j]);
			}
		}
	}

	/* compile the model */
	IloModel model(env);
	for (int i = 0; i < n_vertices; ++i) {
		if (i < n_terminals) {
			for (int j = 0; j < n_terminals; ++j) {
				model.add(terminalRanges[i][j]);
			}
		}
		model.add(vertexRanges[i]);
		for (int j = 0; j < n_vertices; ++j) {
			for (int k = 0; k < n_terminals; ++k) {
				if (this->edge_matrix[i][j] == true) {
					model.add(l1Ranges1[i][j][k]);
					model.add(l1Ranges2[i][j][k]);
				}
			}
		}
	}
	model.add(obj);

	/* solve the model */
	IloCplex solver(model);
	try {
		solver.solve();
	}
	catch (IloException &ex) {
		cerr << ex << endl;
		return -1;
	}
	cout << solver.getObjValue() << endl;
	/* save results*/
	for (int i = 0; i < n_vertices; ++i) {
		cout << "vertex " << i << " : ";
		for (int j = 0; j < n_terminals; ++j) {
			this->simplex_vertices[i][j] = solver.getValue(u[i][j]);
			cout << this->simplex_vertices[i][j] << ' ';
		}
		cout << endl;
	}

	if (CompareDoubleUlps(solver.getObjValue(), 0.0) == 0) {
		return 0.0;
	}
	else
		return solver.getObjValue();
}

double MultiwayCut::post_process(void)
{
	removed_edge = (bool **)malloc(sizeof(bool) * n_vertices * n_vertices);

	for (int i = 0; i < n_vertices; i++) {
		removed_edge[i] = (bool *)malloc(sizeof(bool) * n_vertices);
		for (int j = 0; j < n_vertices; j++) {
			removed_edge[i][j] = false;
		}
	}

	/* For each edge e(u,v), if u and v are not belonged to the same terminal, remove that edge */
	for (int i = 0; i < n_vertices; i++) {
		for (int j = i + 1; j < n_vertices; j++) {
			if (edge_matrix[i][j] == true && assigned_terminal[i] != assigned_terminal[j]) {
				removed_edge[i][j] = true;
			}
		}
	}

	/* For each edge e(u,v), if u and v are not belonged to the same terminal, remove that edge */
	for (int i = 0; i < n_vertices; i++) {
		for (int j = i+1; j < n_vertices; j++) {
			if (edge_matrix[i][j] == true && assigned_terminal[i] != assigned_terminal[j]) {
				removed_edge[i][j] = true;
			}
			else
			{
				removed_edge[i][j] = false;
			}
		}
	}

	/* Calculate sum of the removed edges' weights */
	double sum = 0.0;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = i + 1; j < n_vertices; j++) {
			if (removed_edge[i][j] == true) {
				sum += weight_matrix[i][j];
			}
		}
	}

	//free(assigned_terminal);
	//free(removed_edge);
	return sum;
}

double MultiwayCut::rounding_alg_exp(void)
{
	cout << "Exponential Clock Algorithm" << endl;

	/* Exponential Clock - Terminal sampling */
	double *terminal_clock;
	terminal_clock = new double[n_terminals];
	/* generation of the expoential clocks of the terminals */
	std::default_random_engine generator;

	for (int i = 0; i < n_terminals; ++i) {  // i for terminals
		std::exponential_distribution<double> distribution(1.0);
		//double number = distribution(generator);
		//exponential_clock[i] = y_i*exp(double(-1 * y_i*number));
		terminal_clock[i] = distribution(generator);  // terminal_clock[terminals[i]] denotes ith terminal's exponential clock, terminal[i] denotes ith terminal
	}


	/* Traversing all simplex vertices, find minimum Zi / Ui. then assign the vertice to ith terminal */
	for (int j = 0; j < n_vertices; ++j) {  // j for vertices
		double min = DBL_MAX;
		for (int i = 0; i < n_terminals; ++i) {
			if ((terminal_clock[i] / simplex_vertices[j][i]) < min) {  // Zi / ui
				min = (terminal_clock[terminals[i]] / simplex_vertices[j][i]);  // i = terminal index, terminasl[i] = ith terminal's vertex number

				assigned_terminal[j] = terminals[i];
			}
		}
	}

	//free(terminal_clock);

	return post_process();

}

double MultiwayCut::rounding_alg_dist(void)
{
	cout << "Distortion Algorithm" << endl;

	double r = (double)rand() / RAND_MAX;

	for (int i = 0; i < n_vertices; i++) {
		assigned_terminal[i] = -1;  // -1 is the undefined symbol
	}

	for (int i = 0; i < n_terminals - 1; i++) {
		for (int j = 0; j < n_vertices; j++) {
			if (assigned_terminal[j] == -1 && r < simplex_vertices[j][i] * simplex_vertices[j][i]) {
				assigned_terminal[j] = terminals[i];
			}
		}
	}

	for (int i = 0; i < n_vertices; i++) {
		if (assigned_terminal[i] == -1)
			assigned_terminal[i] = terminals[n_terminals - 1];
	}

	return post_process();
}

double MultiwayCut::rounding_alg(void)
{

	assigned_terminal = (int *)malloc(sizeof(int) * n_vertices);
	double rounded_solution = 0.0;
	double r = (double)rand() / RAND_MAX;
	cout << "r-value : " << r << endl;


	if (r <= 2.0 / 3.0) {
		rounded_solution = MultiwayCut::rounding_alg_exp();
	}
	else {
		rounded_solution = MultiwayCut::rounding_alg_dist();
	}

	return rounded_solution;
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

int MultiwayCut::get_n_vertices() {
	return n_vertices;
}

bool** MultiwayCut::get_removed_edge() {
	return removed_edge;
}