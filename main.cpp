#include "MultiwayCut.h"

int main(void) {

	MultiwayCut a;
	cout << "relaxed_solution: " << a.LP_solver() << endl;
	cout << "rounded_solution: " << a.rounding_alg() << endl;

	return 0;
}