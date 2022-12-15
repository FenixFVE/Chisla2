#pragma once

#include "diag_matrix.hpp"
#include "block.hpp"
#include <cmath>


// jacobi
double step(
	const DiagMatrix& diag_matrix,
	const vector<double>& current_x,
	vector<double>& next_x,
	const vector<double>& F,
	double relaxation) {

	auto& matrix = diag_matrix.diag;
	auto& indexes = diag_matrix.indexes;
	int n = diag_matrix.n;
	double residual = 0.0;
	double F2_sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < 7; j++)
		{
			if (indexes[j] + i >= 0 && indexes[j] + i < n)
			{
				sum += matrix[j][i] * current_x[indexes[j] + i];
			}
		}
		residual += (F[i] - sum) * (F[i] - sum);
		F2_sum += F[i] * F[i];
		next_x[i] = current_x[i] + relaxation / matrix[3][i] * (F[i] - sum);
	}
	return sqrt(residual / F2_sum);
}


// gauss
double step(
	const DiagMatrix& diag_matrix,
	vector<double>& current_x,
	const vector<double>& F,
	double relaxation) {

	auto& matrix = diag_matrix.diag;
	auto& indexes = diag_matrix.indexes;
	int n = diag_matrix.n;
	double residual = 0.0;
	double F2_sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < 7; j++)
		{
			if (indexes[j] + i >= 0 && indexes[j] + i < n)
			{
				sum += matrix[j][i] * current_x[indexes[j] + i];
			}
		}
		residual += (F[i] - sum) * (F[i] - sum);
		F2_sum += F[i] * F[i];
		current_x[i] += relaxation / matrix[3][i] * (F[i] - sum);
	}
	return sqrt(residual / F2_sum);
}

//block
double block_step(
	DiagMatrix& block_diag_matrix, 
	vector<double>& current_x, 
	const vector<double>& F,
	double relaxation,
	int block_size,
	vector<double>& r) {
	
	auto& n = block_diag_matrix.n;
	auto& matrix = block_diag_matrix.diag;
	auto& indexes = block_diag_matrix.indexes;

	double residual = 0.0;
	double sum_of_sq_vec_F = 0.0;
	auto n_blocks = n / block_size;
	for (int i = 0; i < n_blocks; i++) {
		auto k0 = i * block_size;
		auto k1 = (i + 1) * block_size;
		calculate_block_part(block_diag_matrix, current_x, r, k0, k1);
		auto bi = 0;
		for (int j = k0; j < k1; j++, bi++) {
			double sum = 0.0;
			for (int k = 0; k < 7; k++) {
				if (indexes[k] + j >= 0 && indexes[k] + j < n) {
					if (indexes[k] + j < k0 || indexes[k] + j >= k1) {
						sum += matrix[k][j] * current_x[indexes[k] + j];
					}
				}
			}
			r[bi] = F[j] - (r[bi] + sum);
			residual += (r[bi]) * (r[bi]);
			r[bi] *= relaxation;
			sum_of_sq_vec_F += F[j] * F[j];
		}
		solve_SLAE(block_diag_matrix, r, relaxation, k0, k1, block_size);
		bi = 0;
		for (int j = k0; j < k1; j++, bi++)
		{
			current_x[j] += r[bi];
		}
	}
	return sqrt(residual / sum_of_sq_vec_F);
}