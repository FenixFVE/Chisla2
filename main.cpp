
#include "diag_matrix.hpp"
#include "gauss.hpp"
#include "jacobi.hpp"
#include "block.hpp"
#include "file_manager.hpp"

using namespace std;

int main() {
	try {
		int size, block_size, max_steps;
		double relaxation, eps;

		read_parameters(size, relaxation, eps, max_steps, block_size, "parameters.txt");
		cout << "Relaxation: " << relaxation << " Eps: " << eps << " Max steps: " << max_steps << " Block_size: " << block_size << "\n\n";

		vector<double> temporal(size, 0.0);
		auto x_jacobi = read_vector("initial_vector.txt");
		auto x_gauss = x_jacobi;
		auto x_block = x_jacobi;
		auto F = read_vector("vector_Fa.txt");
		DiagMatrix matrix("matrix_A.txt");

		jacobi(matrix, F, x_jacobi, relaxation, eps, max_steps, temporal);
		write_vector("jacobi_output.txt", x_jacobi);
		cout << '\n';

		gauss(matrix, F, x_gauss, relaxation, eps, max_steps);
		write_vector("gauss_output.txt", x_gauss);
		cout << '\n';

		block(matrix, F, x_block, relaxation, eps, max_steps, block_size, temporal);
		write_vector("block_output.txt", x_block);
	}
	catch (exception exc) {
		cout << "Error: " << exc.what() << '\n';
	} 
	return 0;
}
