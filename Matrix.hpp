
#ifndef EX3_H

#include <iostream>
#include <vector>
#include "Complex.h"
#include <algorithm>
#include <unordered_map>
#include <future>
#include <assert.h>
#include <functional>
#include <assert.h>


#define INVALID_MATRIX_SIZE "Matrix size isn’t equal, getRows or columns"
#define MATRICES_MUST_MATCH "Matrices's size isn’t must match!"
#define BAD_MULTIPLICATION "col of Matrix1 and row of the Matrix2 are not equal"
#define OUT_OF_RANGE  "The requested index is not inside the matrix boundaries"
#define PARALLEL_MSG "Generic Matrix mode changed to parallel mode.\n"
#define NON_PARALLEL_MSG "Generic Matrix mode changed to non-parallel mode.\n"

using namespace std;

template<typename T>
/**
 *
 * @tparam T
 */
class Matrix
{
    size_t row;
    size_t col;
    vector<T> matrix;
    bool static isParallel;

public:

    using const_iterator = typename vector<T>::const_iterator ;

// -----------------------------------Declarations------------------------------------

    //----------- Constructors-------------------------------------


    /**
     * Default constructor
     * @return a 1x1 matrix, with 0 at the cell.
     */
    Matrix();

    /**
     * Constructor of a matrix with a given rows and columns.
     * @param rows
     * @param cols
     */
    Matrix(size_t rows, size_t cols);

    /**
     * Copy constructor
     * @param object
     */
    Matrix(const Matrix &object);

    /**
     *Move constructor- optional
     * @param otherMatrix
     */
    Matrix(Matrix && otherMatrix);

    /**
     * Constractor of a matrix with a given size and values.
     * @param rows
     * @param cols
     * @param cells
     */
    Matrix(size_t rows, size_t cols, const vector<T> &cells);

    /**
     * Destructor
     */
    ~Matrix()
    {};

    //-----------------Operators-----------------------------------------


    /**
     *
     * @param otherMatrix
     * @return
     */
    Matrix &operator=(const Matrix<T> &otherMatrix);

    /**
     *
     * @return
     */
    size_t size() const;

    /**
     *
     * @param otherMatrix
     * @return
     */
    Matrix<T> operator+(const Matrix<T> &otherMatrix) const;

    /**
     *
     * @param otherMatrix
     * @return
     */
    Matrix<T> operator-(const Matrix<T> &otherMatrix) const;

    /**
     *
     * @param otherMatrix
     * @return
     */
    Matrix<T> operator*(const Matrix<T> &otherMatrix) const;

    /**
     *
     * @param otherMatrix
     * @return
     */
    bool operator!=(const Matrix<T> &otherMatrix) const;

    /**
     *
     * @param otherMatrix
     * @return
     */
    bool operator==(const Matrix<T> &otherMatrix) const;

    /**
     *
     * @param _row
     * @param _col
     * @return
     */
    T operator()(size_t _row, size_t _col) const ;

    /**
     *
     * @param _row
     * @param _col
     * @return
     */
    T &operator()(size_t _row, size_t _col);

    /**
     *
     * @return
     */
    Matrix<T> trans();

    /**
     *
     * @return
     */
    bool isSquareMatrix();

    /**
     *
     * @return
     */
    typename vector<T>::const_iterator begin() const ;

    /**
     *
     * @return
     */
    typename vector<T>::const_iterator end() const ;

    /**
     *
     * @return
     */
    size_t getRows() const;

    /**
     *
     * @return
     */
    size_t getCols() const;

    /**
     *
     * @param state
     */
    static void setParallel(bool state);

    /**
     *
     * @tparam U
     * @param os
     * @param toPrint
     * @return
     */
    template<typename U>
    friend std::ostream &operator<<(ostream &os, const Matrix<U> &toPrint);

};



    // -----------------------------------Implementations------------------------------------

//-----------------Constructors-----------------------------------------

    template<typename T>
    Matrix<T>::Matrix(): row(1), col(1), matrix(1, T(0)) {}


    template<typename T>
    Matrix<T>::Matrix(const Matrix &object): row(object.row), col(object.col)
    {
        matrix = object.matrix;
    }


    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols): row(rows), col(cols)
    {
        matrix.assign((row * col), T(0));
    }

    template<typename T>
    Matrix<T>::Matrix(Matrix && otherMatrix): row(otherMatrix.row), col(otherMatrix.col)
    {
        this->matrix.swap(otherMatrix.matrix);
        otherMatrix.matrix.clear();
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, const vector<T> &cells): row(rows), col(cols),
                                                                         matrix(cells) {}



//-----------------Operators-----------------------------------------


template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &otherMatrix)
{
    if (&otherMatrix != this)
    {
        row = otherMatrix.row;
        col = otherMatrix.col;
        matrix.clear();
        matrix = otherMatrix.matrix;
    }
    return *this;
}


template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &otherMatrix) const
{
    if (otherMatrix.row != row || otherMatrix.col != col)
    {
        throw out_of_range(INVALID_MATRIX_SIZE);
    }
    Matrix<T> resultedMatrix = Matrix(row, col);
    if (!isParallel)
    {
        for (unsigned int i = 0; i < resultedMatrix.size() ; ++i)
        {
            resultedMatrix[i] = matrix[i] + otherMatrix[i];
        }
    }
    else
    {
        vector<future<void>> parallelVec;
        for (unsigned int i = 0; i < row; ++i)
        {
            parallelVec.emplace_back(async([i, this, &otherMatrix, &resultedMatrix]()
                                           {
                                               for (unsigned int j = 0; j < col; ++j)
                                               {
                                                   resultedMatrix(i, j) = (*this)(i, j) +
                                                           otherMatrix(i, j);
                                               }
                                           }));
            for (auto &k : parallelVec)
            {
                k.wait();
            }
        }
    }
    return resultedMatrix;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &otherMatrix) const
{
    if (otherMatrix.row != row || otherMatrix.col != col)
    {
        throw out_of_range(MATRICES_MUST_MATCH);
    }
    Matrix resultedMatrix = Matrix(row, col);
    for (int i = 0; i < otherMatrix.size(); ++i)
    {
        resultedMatrix.matrix[i] = matrix[i] - otherMatrix.matrix[i];
    }
    return resultedMatrix;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &otherMatrix) const
{
    if (col != otherMatrix.row)
    {
        throw invalid_argument(BAD_MULTIPLICATION);
    }
    Matrix<T> resultedMatrix = Matrix(this->row, otherMatrix.col);
    if (!isParallel)
    {
        for (unsigned int i = 0; i < this->row; ++i)
        {
            for (unsigned int j = 0; j < otherMatrix.col; ++j)
            {
                for (unsigned int k = 0; k < this->col; ++k)
                {
                    resultedMatrix(i, j) += (*this)(i, k) * otherMatrix(k, j);
                }
            }

        }
        return resultedMatrix;
    }

    vector<future<void>> parallelVec;
    for (unsigned int i = 0; i < this->row; ++i)
    {
        parallelVec.emplace_back(async([i, this, &otherMatrix, &resultedMatrix]()
                                       {
                                           for (unsigned int j = 0; j < otherMatrix.col; ++j)
                                           {
                                               for (unsigned int k = 0; k < this->col; ++k)
                                               {
                                                   resultedMatrix(i, j) += ((*this)(i, k) *
                                                           otherMatrix(k, j));
                                               }
                                           }
                                       }));
    }
    for (unsigned int l = 0 ; l < parallelVec.size(); ++l)
    {
        parallelVec[l].wait();
    }
    return resultedMatrix;
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix<T> &otherMatrix) const
{
    if (otherMatrix.row != row || otherMatrix.col != col)
    {
        throw invalid_argument(INVALID_MATRIX_SIZE);
    }
    return (this != otherMatrix);
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &otherMatrix) const
{
    if (otherMatrix.row != row || otherMatrix.col != col)
    {
        throw invalid_argument(INVALID_MATRIX_SIZE);
    }
    if ((row != otherMatrix.row) || (col != otherMatrix.col))
    {
        return false;
    }

    return !(matrix != otherMatrix.matrix);
}

template<typename T>
T Matrix<T>::operator()(size_t _row, size_t _col) const
{
    if (_row > row || _col > col)
    {
        throw out_of_range(OUT_OF_RANGE);
    }

    return matrix[(col * _row) + _col];
}

template<typename T>
T &Matrix<T>::operator()(size_t _row, size_t _col)
{
    if (_row > row || _col > col)
    {
        throw out_of_range(OUT_OF_RANGE);
    }
    return matrix[(col * row) + _col];
}


template<typename T>
size_t Matrix<T>::size() const
{
    return matrix.size();
}

template<typename T>
Matrix<T> Matrix<T>::trans()
{
    Matrix transposedMat = Matrix(col, row);
    for (unsigned int i = 0; i < row; ++i)
    {
        for (unsigned int j = 0; j < col; ++j)
        {
            transposedMat(j, i) = (*this)(i, j);
        }
    }
    return transposedMat;
}

template<typename T>
bool Matrix<T>::isSquareMatrix()
{
    return (row == col);
}

template<typename T>
typename vector<T>::const_iterator Matrix<T>::begin() const
{
    return matrix.begin();
}

template<typename T>
typename vector<T>::const_iterator Matrix<T>::end() const
{
    return matrix.end();
}

template<typename T>
size_t Matrix<T>::getRows() const
{
    return row;
}

template<typename T>
size_t Matrix<T>::getCols() const
{
    return col;
}


template <typename T>
ostream &operator<<(ostream &os, const Matrix<T> &toPrint)
{
    for (int i = 0; i < toPrint.size() ; ++i)
    {
        ((i + 1) % toPrint.col == 0) ? os << toPrint.matrix[i] << '\t' << '\n' : os << toPrint
        .matrix[i] << '\t';
    }
    return os;
}

template <>
inline Matrix<Complex> Matrix<Complex>::trans()
{
    Matrix transCompMat = Matrix(col, row);
    for (unsigned int i = 0; i < row; ++i)
    {
        for (unsigned int j = 0; j < col; ++j)
        {
            Complex ans = (*this)(i, j).conj();
            transCompMat(j, i) = ans;
        }
    }
    return transCompMat;
}

template<typename T>
void Matrix<T>::setParallel(bool state)
{
    {
        state ? std::cout << PARALLEL_MSG << endl :
        cout << NON_PARALLEL_MSG << std::endl;
        isParallel = state;
    }
}
#endif // EX3_H