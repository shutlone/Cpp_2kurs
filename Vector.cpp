#include <cmath>
#include "Vector.h"
#include "Matrix.h"

using namespace mat_vec;

Vector::Vector(size_t size, double value): length(size),data(new double[size]) {

    for(size_t i = 0; i < size; i++)
        data[i] = value;
}


Vector::Vector(const mat_vec::Vector &src)
{
    this->length = src.length;
    this->data = new double[length];
    for(size_t i = 0; i < length; i++)
    {
        this->data[i] = src.data[i];
    }
}

Vector& Vector::operator=(const mat_vec::Vector &rhs){
    delete[] this->data; // освобождаем память, чтобы не было утечки памяти
    this->length = rhs.length;
    this->data = new double[length];
    for(size_t i = 0; i < length; i++)
    {
        this->data[i] = rhs.data[i];
    }
    return *this; // возвращает указатель на этот вектор
}

Vector::~Vector()
{
    delete[] data;
    data = nullptr;
}

size_t Vector::size() const {
    return this->length;
}

double Vector::operator[](size_t n) const {
    return this->data[n];
}

double& Vector::operator[](size_t n) {
  return this->data[n];
}

double Vector::norm() const {
    Vector V(length);
  for(size_t i = 0; i < length; i++)
  {
    V.data[i] = this->data[i];
  }
    double L_norm = 0;
    for(size_t i = 0; i < V.length; i++)
       L_norm += V[i] * V[i];
    L_norm = sqrt(L_norm);
    return L_norm;
}

Vector Vector::normalized() const {
  Vector V(length);
  double length_vec = 0;

  for(size_t i = 0; i < length; i++)
  {
    V.data[i] = this->data[i];
  }

  for(size_t i = 0; i < length; i++)
    length_vec += V.data[i] * V.data[i];

  length_vec = sqrt(length_vec);

  for(size_t i = 0; i < length; i++)
    V.data[i] = V.data[i] / length_vec;

  return V;
}

void Vector::normalize() {
  double length_vec = 0;
  for(size_t i = 0; i < length; i++)
    length_vec += data[i] * data[i];

  length_vec = sqrt(length_vec);

  for(size_t i = 0; i < length; i++)
    data[i] = data[i] / length_vec;
}

Vector Vector::operator+(const mat_vec::Vector &rhs) const {
  Vector V(length);
    for(size_t i = 0; i < length; i++){
      V.data[i] = this->data[i] + rhs.data[i];
    }
  return V;
}

Vector& Vector::operator+=(const mat_vec::Vector &rhs) {
  return *this = *this + rhs;
}

Vector Vector::operator-(const mat_vec::Vector &rhs) const {
  Vector V(length);
  for(size_t i = 0; i < length; i++){
    V.data[i] = this->data[i] - rhs.data[i];
  }
  return V;
}

Vector& Vector::operator-=(const mat_vec::Vector &rhs) {
  return *this = *this - rhs;
}

Vector Vector::operator^(const mat_vec::Vector &rhs) const {
  Vector V(length);
  for(size_t i = 0; i < length; i++)
    V[i] = this->data[i] * rhs.data[i];
  return V;
}

Vector& Vector::operator^=(const mat_vec::Vector &rhs) {
  return *this = *this ^ rhs;
}

double Vector::operator*(const mat_vec::Vector &rhs) const {
  double scalar = 0;
  for(size_t i = 0; i < length; i++)
    scalar += this->data[i] * rhs.data[i];
  return scalar;
}

Vector Vector::operator*(double k) const {
  Vector V(length);
  for(size_t i = 0; i < length; i++)
    V.data[i] = data[i] * k;
  return V;
}

Vector& Vector::operator*=(double k) {
  for(size_t i = 0; i < length; i++)
    data[i] = data[i] * k;
  return *this;
}

Vector Vector::operator/(double k) const {
  Vector V(length);
  for(size_t i = 0; i < length; i++)
    V.data[i] = data[i] / k;
  return V;
}

Vector& Vector::operator/=(double k) {
  for(size_t i = 0; i < length; i++)
    data[i] = data[i] / k;
  return *this;
}

Vector Vector::operator*(const mat_vec::Matrix &mat) const {
  size_t rows = mat.shape().first;
  size_t cols = mat.shape().second;
  double sum = 0;
  Vector V(cols);

  for(size_t m = 0; m < cols; m++) {
    for (size_t n = 0; n < rows; n++)
      sum += this->data[n] * mat.get(m, n);
    V.data[m] = sum;
    sum = 0;
  }

  return V;
}

Vector& Vector::operator*=(const mat_vec::Matrix &mat) {
  return *this = *this * mat;
}

bool Vector::operator==(const mat_vec::Vector &rhs) const {
  if(this->length != rhs.length)
    return false;
  if(this->length == rhs.length)
    for(size_t i = 0; i < length; i++){
      if(data[i] != rhs.data[i])
        return false;
    }
  return true;
}

bool Vector::operator!=(const mat_vec::Vector &rhs) const {
  if(this->length != rhs.length)
    return true;
  if(this->length == rhs.length)
    for(size_t i = 0; i < length; i++){
      if(this->data[i] != rhs.data[i])
        return true;
    }
  return false;
}

Vector mat_vec::operator*(double k, const Vector &v){
  return v * k;
};
