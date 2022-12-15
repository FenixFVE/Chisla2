#pragma once

#include "file_manager.hpp"

using namespace std;

class DiagMatrix {
public:
	int n;
	int m;
	vector<int> indexes;
	vector<vector<double>> diag;
	DiagMatrix(const string& matrix_file) {
		try {
			auto file = ifstream(matrix_file);
			if (!file.good()) throw 0;
			file >> n >> m;

			indexes = vector<int>(7);
			for (int i = 0; i < 7; ++i) {
				file >> indexes[i];
			}

			diag = vector<vector<double>>(7, vector<double>(n));
			for (int i = 0; i < 7; ++i) {
				for (int j = 0; j < n; ++j) {
					file >> diag[i][j];
				}
			}
		}
		catch (...) {
			throw exception(format("Fail to read matrix: \"{}\"", matrix_file).c_str());
		}
	}

	void LU_decompositon(int block_size) {
		int N = n / block_size;
		for (int i = 0; i < N; ++i) {
			int k0 = i * block_size;
			int k1 = (i + 1) * block_size;
			for (int j = k0 + 1; j < k1; ++j) {
				diag[4][j - 1] /= diag[3][j - 1];
				diag[3][j] -= diag[4][j - 1] * diag[2][j];
			}
		}
	}
};


double cond(const vector<double>& x, double residual)
{
	int n = x.size();
	double x_star_norm = 0.0, x_minus_x_star_norm = 0.0;
	for (int i = 0; i < n; ++i) {
		double i1 = static_cast<double>(i + 1);
		x_star_norm += i1 * i1;
		x_minus_x_star_norm += (x[i] - i1) * (x[i] - i1);
	}
	x_star_norm = sqrt(x_star_norm);
	x_minus_x_star_norm = sqrt(x_minus_x_star_norm);
	auto error = x_minus_x_star_norm / x_star_norm;
	return error / residual;
}