#include "MultiwayCut.h"

int main(void) {
	double LP;
	double RS;
	
	/* loop until different between relaxed solution & rounded solution
	do {
		MultiwayCut a;
		LP = a.LP_solver();
		RS = a.rounding_alg();
		cout << "relaxed_solution: " << LP << endl;
		//cout << "optimal_solution: " << a.get_optimal_solution() << endl;
		cout << "rounded_solution: " << RS << endl;
	} while (CompareDoubleUlps(LP,RS) == 0);
	*/


	MultiwayCut a;
	LP = a.LP_solver();
	RS = a.rounding_alg();
	cout << "relaxed solution: " << LP << endl;
	cout << "rounded solution: " << RS << endl;


	cout << "Different!" << endl;
	return 0;
}