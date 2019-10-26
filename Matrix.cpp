#include "Matrix.h"
#include "Vector.h"
#include <queue>
#include <cmath>
#include <iostream>

using namespace mat_vec;

Matrix::Matrix(size_t size, double value): rows(size),cols(size),data(new double*[size]) {

    for (size_t i = 0; i < size; i++) {
      data[i] = new double[size];
      for(size_t j = 0; j < size; j++)
        data[i][j] = value;
    }

}


Matrix Matrix::eye(size_t size) {
  Matrix M(size, size, 0);
  for(size_t i = 0; i < M.cols; i++)
    M.data[i][i] = 1;
  return M;
}

Matrix::Matrix(size_t rows, size_t cols, double value): rows(rows),cols(cols),data(new double *[rows]) {

    for(size_t i = 0; i < rows; i++) {
      data[i] = new double[cols];
      for (size_t j = 0; j < cols; j++)
        data[i][j] = value;
    }

}

Matrix::Matrix(const mat_vec::Matrix &src) {
  this->rows = src.rows;
  this->cols = src.cols;
  this->value = src.value;
  if(src.data != nullptr) {

    this->data = new double *[rows];

    for (size_t i = 0; i < rows; i++) {
      data[i] = new double[cols];
      for (size_t j = 0; j < cols; j++)
        data[i][j] = src.data[i][j];
    }
  } else
    this->data = nullptr;

}

Matrix& Matrix::operator=(const mat_vec::Matrix &rhs) {
  delete[] this->data;
  this->rows = rhs.rows;
  this->cols = rhs.cols;

  data = new double*[rows];
  for(size_t i = 0; i < rows; i++)
    data[i] = new double[cols];

  for(size_t i = 0; i < rows; i++)
    for(size_t j = 0; j < cols; j++)
      data[i][j] =rhs.data[i][j];
  return *this;
}

Matrix::~Matrix() {
  for(size_t i = 0; i < rows; i++)
    delete []data[i];

  delete []data;

}

void Matrix::reshape(size_t rows, size_t cols) {
  std::queue<double> q;

  for(size_t i = 0; i < this->rows; i++)
    for(size_t j = 0; j < this-> cols; j++)
      q.push(this->data[i][j]);

  for(size_t i = 0; i < rows; i++)
    delete[] data[i];

  delete[] data;
  data = nullptr;

  this->rows = rows;
  this->cols = cols;
  data = new double *[rows];
  for (size_t i = 0; i < rows; i++)
    data[i] = new double[cols];

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      this->data[i][j] = q.front();


}

std::pair<size_t, size_t> Matrix::shape() const {
  std::pair<size_t ,size_t > temp;
  temp.first = this->rows;
  temp.second = this->cols;

  return temp;
}

double Matrix::get(size_t row, size_t col) const {
  return data[row][col];
}

Matrix Matrix::operator+(const mat_vec::Matrix &rhs) const {
  Matrix M(rows,cols);
  for(size_t i = 0; i < rows; i++){
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = this->data[i][j] + rhs.data[i][j];
  }
  return M;
}

Matrix& Matrix::operator+=(const mat_vec::Matrix &rhs) {
  return *this = *this + rhs;
}

Matrix Matrix::operator-(const mat_vec::Matrix &rhs) const {
  Matrix M(rows,cols);
  for(size_t i = 0; i < rows; i++){
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = this->data[i][j] - rhs.data[i][j];
  }

  return M;
}

Matrix& Matrix::operator-=(const mat_vec::Matrix &rhs) {
  return *this = *this - rhs;
}

Matrix Matrix::operator*(const mat_vec::Matrix &rhs) const {
  Matrix M(this->rows, rhs.cols);
  for(size_t i = 0; i < this->rows; i++) {
    for (size_t j = 0; j < rhs.cols; j++) {
      for (size_t k = 0; k < this->cols; k++)
        M.data[i][j] = (M.data[i][j] + this->data[i][k] * rhs.data[k][j]);
    }
  }
  return M;
}

Matrix& Matrix::operator*=(const mat_vec::Matrix &rhs) {
  return *this = *this * rhs;
}

Matrix Matrix::operator*(double k) const {
  Matrix M(rows,cols,0.0);
  for(size_t i = 0; i < rows; i++)
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = M.data[i][j] * k;

  return M;
}

Matrix& Matrix::operator*=(double k) {
  return *this = *this * k;
}

Matrix Matrix::operator/(double k) const {
  Matrix M(rows,cols);
  for(size_t i = 0; i < rows; i++)
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = M.data[i][j] / k;

  return M;
}

Matrix& Matrix::operator/=(double k) {
  return *this = *this / k;
}

Matrix Matrix::transposed() const {
  Matrix M(this->rows,this->cols);
  for(size_t i = 0; i < rows; i++)
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = this->data[i][j];
  M.transpose();
  return M;
}

void Matrix::transpose() {
  double **M = new double*[cols];
  for(size_t i = 0; i < cols; i++){
    M[i] = new double[rows];
    for(size_t j = 0; j < rows; j++){
      M[i][j] = data[j][i];
    }
  }

  for(size_t i = 0; i < rows; i++)
    delete []data[i];

  delete []data;

  double temp = rows; rows = cols; cols = temp;
  data = new double*[rows];
  for(size_t i = 0; i < rows; i++)
  {
    data[i] = new double[cols];
    for(size_t j = 0; j < cols; j++)
    {
      data[i][j] = M[i][j];
    }
  }
  for(size_t i = 0; i < rows; i++)
  {
    delete []M[i];
  }
  delete []M;
}

double Matrix::det() const {
  Matrix M(this->rows,this->cols);

  for(size_t i = 0; i < rows; i++)
    for(size_t j = 0; j < cols; j++)
      M.data[i][j] = this->data[i][j];

  const double EPS = 1E-9;
  double determinant = 1;
    for (size_t i = 0; i < M.rows; i++) {
      size_t k = i;
      for (size_t j = i + 1; j < M.cols; j++)
        if (std::abs(M.data[j][i]) > std::abs(M.data[k][i]))
          k = j;
      if (std::abs(M.data[k][i]) < EPS) {
        return 0;
      }
      for (size_t j = 0; j < M.cols; j++) {
        double temp = M.data[i][j];
        M.data[i][j] = M.data[k][j];
        M.data[k][j] = temp;
      }
      if (i != k)
        determinant = -determinant;
      determinant = determinant * M.data[i][i];
      for (size_t j = i + 1; j < M.rows; j++)
        M.data[i][j] = M.data[i][j] / M.data[i][i];
      for (size_t j = 0; j < M.rows; j++)
        if (j != i && std::abs(M.data[j][i]) > EPS)
          for (size_t m = i + 1; m < M.rows; k++)
            M.data[j][m] = M.data[j][m] - M.data[i][m] * M.data[j][i];

    }
    return determinant;


}

Matrix Matrix::inv() const{
  Matrix M(rows, cols, 0.0);
  for(size_t i = 0; i < rows; i++){
    for(size_t j = 0; j < rows; j++){
      M.data[i][j] = data[i][j];
    }
  }
  double det = M.det();

  for (size_t i = 0; i< rows; i++){
    if (data[i][i] == 0)
      for (size_t j = i+1; j< rows; j++){
        if (data[j][i]==1){
          for (size_t k =0; k< 2*rows; k++){
            double c = data[j][k];
            data[j][k] = data[i][k];
            data[i][k] = c;
          }
          break;
        }
      }
    for (size_t k = i+1; k < rows; k++){
      if (data[k][i] == 1){
        for (size_t j = 0; j < 2*rows; j++){
          data[k][j] *= data[i][j];
        }
      }
    }
  }
  for (size_t i = rows-1; i >= 0; i--){
    for (size_t k = i-1; k >= 0; k--){
      if (data[k][i] == 1){
        for (size_t j = 0; j < 2*rows; j++){
          data[k][j] *= data[i][j];
        }
      }
    }
  }

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0, k = rows; j < rows; j++, k++)
      M.data[i][j] = data[i][k];

  Matrix E = M.eye(rows);

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < rows; j++)
      for (size_t k = 0; k < rows; k++)
        E.data[i][j] *= data[i][k] * data[k][j];

  return M;
}

Vector Matrix::operator*(const mat_vec::Vector &vec) const {
  Vector V(rows);
  double sum = 0;
  for(size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++)
      sum += this->data[i][j] * vec[j];
    V[i] = sum;
    sum = 0;
  }

  return V;
}

bool Matrix::operator==(const mat_vec::Matrix &rhs) const {
  if(this->rows != rhs.rows || this->cols != rhs.cols)
    return false;
  for(size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++)
      if (data[i][j] != rhs.data[i][j])
        return false;
  }

  return true;
}

bool Matrix::operator!=(const mat_vec::Matrix &rhs) const {
  if(this->rows != rhs.rows || this->cols != rhs.cols)
    return true;
  for(size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++)
      if (data[i][j] != rhs.data[i][j])
        return true;
  }

  return false;
}
