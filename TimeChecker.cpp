//
// Created by nadav.meidan on 1/19/19.
//

#include <iostream>
#include "Matrix.hpp"
#include <eigen3/Eigen/Dense>
#include <stack>
#include <ctime>
#include <chrono>

using namespace std;
using namespace Eigen;

const int ARGS_NUM = 2;

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[])
{
    if (argc < ARGS_NUM)
    {
        cerr << "Not enough arguments";
        return EXIT_FAILURE;
    }

    stack<chrono::time_point<chrono::system_clock>> tictoc_stack;

    unsigned int size = (unsigned int)stoi(argv[1]);
    cout << "size " << size << endl;

    // initializing eigen matrices
    MatrixXd eigenMat1 = MatrixXd::Random(size, size);
    MatrixXd eigenMat2 = MatrixXd::Random(size, size);

    // eigen time measure for +\* matrices operations
    tictoc_stack.push(chrono::system_clock::now());
    eigenMat1 * eigenMat2;
    chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - tictoc_stack.top();
    cout << "eigen mult " << elapsed_seconds.count() << endl;
    tictoc_stack.pop();

    tictoc_stack.push(chrono::system_clock::now());
    eigenMat1 + eigenMat2;
    elapsed_seconds = chrono::system_clock::now() - tictoc_stack.top();
    cout << "eigen add " << elapsed_seconds.count() << endl;
    tictoc_stack.pop();

    // initializing Matrix<int> matrices
    vector<int> vector(size * size, 1);
    Matrix<int> matrix(size, size, vector);

    // Matrix<int> time measure for +\* matrices operations
    tictoc_stack.push(chrono::system_clock::now());
    matrix * matrix;
    elapsed_seconds = chrono::system_clock::now() - tictoc_stack.top();
    cout << "matlib mult " << elapsed_seconds.count() << endl;
    tictoc_stack.pop();

    tictoc_stack.push(chrono::system_clock::now());
    matrix + matrix;
    elapsed_seconds = std::chrono::system_clock::now() - tictoc_stack.top();
    cout << "matlib add " << elapsed_seconds.count() << endl;
    tictoc_stack.pop();

    return 0;
}