#pragma once

#include "file_manager.hpp"
#include "diag_matrix.hpp"
#include "step.hpp"

double block_step(DiagMatrix&,vector<double>&,const vector<double>&,double,int,vector<double>&);

void block(DiagMatrix& block_diag_matrix, 
    const vector<double>& F, 
    vector<double>& x, 
    double relaxation, 
    double eps, 
    int max_steps,
    int block_size,
    vector<double>& r) {

    cout << "Block Relaxation method start" << endl;
    block_diag_matrix.LU_decompositon(block_size);
    auto residual = DBL_MAX;
    for (int i = 1; i <= max_steps && residual > eps; i++)
    {
        residual = block_step(block_diag_matrix, x, F, relaxation, block_size, r);
        cout << "Iteration number: " << i << " Residual: " << fixed << setprecision(16) << residual << "\r";
    }
    auto con = cond(x, residual);
    cout << "\nCond: " << fixed << setprecision(16) << con << "\nBlock Relaxation method end\n";
}

void calculate_block_part(
    DiagMatrix& block_diag_matrix, 
    vector<double>& x, 
    vector<double>& r, 
    int k0, int k1) {

    auto& n = block_diag_matrix.n;
    auto& matrix = block_diag_matrix.diag;
    auto& indexes = block_diag_matrix.indexes;
    int k = 0;
    for (int i = k0; i < k1; i++, k++) {
        double sum = 0.0;
        for (int j = 3; j < 5; j++) {
            if (indexes[j] + i >= 0 && indexes[j] + i < n) {
                if (indexes[j] + i >= k0 && indexes[j] + i < k1) {
                    if (j == 3) {
                        sum += x[indexes[j] + i];
                    } else {
                        sum += matrix[j][i] * x[indexes[j] + i];
                    }
                }
            }
        }
        r[k] = sum;
    }
    auto buf = r;
    k = 0;
    for (int i = k0; i < k1; i++, k++) {
        double sum = 0.0;
        for (int j = 2; j < 4; j++) {
            if (indexes[j] + i >= 0 && indexes[j] + i < n) {
                if (indexes[j] + i >= k0 && indexes[j] + i < k1) {
                    sum += matrix[j][i] * buf[indexes[j] + k];
                }
            }
        }
        r[k] = sum;
    }
}

void solve_SLAE(
    DiagMatrix& block_diag_matrix, 
    vector<double>& y, 
    double relaxation, 
    int k0, int k1, int block_size) {
    
    auto& matrix = block_diag_matrix.diag;

    int j = 0;
    y[j] = y[j] / matrix[3][k0]; j++;
    for (int i = k0 + 1; i < k1; i++, j++) {
        y[j] = (y[j] - matrix[2][i] * y[j - 1]) / matrix[3][i];
    }

    j = block_size - 2;
    for (int i = k1 - 2; i >= k0; i--, j--) {
        y[j] -= matrix[4][i] * y[j + 1];
    }
}