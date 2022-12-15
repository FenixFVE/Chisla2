#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <exception>
#include <format>
#include <iomanip>
#include "diag_matrix.hpp"

using namespace std;

void read_parameters(int& size, double& relaxation, double& eps, int& max_steps, int& block_size, const string& parameters_file) {
	try {
		auto file = ifstream(parameters_file);
		if (!file.good()) throw 0;
		file >> size >> relaxation >> eps >> max_steps >> block_size;
	}
	catch(...) {
		throw exception(format("Fail to read parameters: \"{}\"", parameters_file).c_str());
	}
}

vector<double> read_vector(const string& vector_file) {
	try {
		int size;
		auto file = ifstream(vector_file);
		if (!file.good()) throw 0;
		file >> size;
		vector<double> vec(size);
		for (int i = 0; i < size; ++i) {
			file >> vec[i];
		}
		return vec;
	}
	catch (...) {
		throw exception(format("Fail to read vector: \"{}\"", vector_file).c_str());
	}
}

void write_vector(const string& vector_file, const vector<double>& vec) {
	try {
		auto file = ofstream(vector_file);
		if (!file.good()) throw 0;
		file << fixed << setprecision(16);
		for (int i = 0; i < vec.size(); ++i) {
			file << vec[i] << '\n';
		}
	} 
	catch (...) {
		throw exception(format("Fail to write vector: \"{}\"", vector_file).c_str());
	}
}
