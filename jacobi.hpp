#pragma once

#include "file_manager.hpp"
#include "diag_matrix.hpp"
#include "step.hpp"

using namespace std;

void jacobi(
	const DiagMatrix& diag_matrix,
	const vector<double>& F,
	vector<double>& x,
	double relaxation,
	double eps,
	int max_steps,
	vector<double>& new_x) {

	cout << "Jacobi method start\n";

	double residual = DBL_MAX;
	for (int i = 1; i <= max_steps && residual > eps; ++i) {
		residual = step(diag_matrix, x, new_x, F, relaxation);
		x = new_x;
		cout << "Iteration: " << i << " Residual: " << fixed << setprecision(16) << residual << '\r';
	}
	double c = cond(x, residual);
	cout << "\nCond: " << fixed << setprecision(16) << c << "\nJacobi method end\n";
}