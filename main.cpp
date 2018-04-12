#include "MultiwayCut.h"

int main(void) {
	double LP;
	double RS;
	int iteration = 0;
	/* loop until different between relaxed solution & rounded solution	*/
	do {
		srand((unsigned)time(NULL) + (unsigned)iteration * 10);
		cout << ++iteration << "th case" << endl << endl;
		MultiwayCut *a = new MultiwayCut;
		LP = a->LP_solver();
		RS = a->rounding_alg();
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