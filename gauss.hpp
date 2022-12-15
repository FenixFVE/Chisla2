#pragma once

#include "diag_matrix.hpp"
#include "file_manager.hpp"
#include "step.hpp"

void gauss(
    const DiagMatrix& diag_matrix, 
    const vector<double>& F,
    vector<double>& x, 
    double relaxation, 
    double eps, 
    int max_steps) {
    cout << "Gauss-Seidel method start" << endl;

    auto residual = DBL_MAX;
    for (int i = 1; i <= max_steps && residual > eps; i++)
    {
        residual = step(diag_matrix, x, F, relaxation);
        cout << "Iteration: " << i << " Residual: " << fixed << setprecision(16) << residual << "\r";
    }
    auto c = cond(x, residual);
    cout << endl << "Cond: " << fixed << setprecision(16) << c << "\nGauss-Seidel method end\n";
}