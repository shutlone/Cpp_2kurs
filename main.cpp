#define CATCH_CONFIG_MAIN

#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Base.h"
#include "catch.hpp"



using namespace mat_vec;

TEST_CASE("Vector"){
  Vector rv(3, 2);
  Vector cv(3, 3.0);
  REQUIRE(rv[0] * cv[0] == 6.0);

  Vector cc(rv);
  REQUIRE(cc[2] == 2.0);

  cc = cv;
  REQUIRE(cc[2] == 3.0);

  Vector b1(4,2);
  Vector b2(4,10);
  REQUIRE(5.0 * b1 == b2);

  Vector c1(2);
  c1[0] = 3;
  c1[1] = 4;
  Vector c2 = c1.normalized();
  REQUIRE(c2[0] == 0.6);
  REQUIRE(c2[1] == 0.8);

  REQUIRE(c2 != c1);

  Vector p1(6,2);
  Vector p2(6,7);
  Vector p3(6,9);

  REQUIRE(p1 + p2 == p3);

  p1 += p2;

  REQUIRE(p1 == p3);

  Vector m1(6,2);
  Vector m2(6,7);
  Vector m3(6,14);

  REQUIRE(m1 * m2 == 84);

  m1 ^= m2;

  REQUIRE(m1 == m3);

  Vector s1(7,3);
  Vector s2(7,4);
  Vector s3(7,-1);

  REQUIRE(s1 - s2 == s3);

  s1 -= s2;

  REQUIRE(s1 == s3);

  Vector ms1(6,5);
  Vector ms2(6,30);

  REQUIRE(ms1 * 6 == ms2);

  ms1 *= 6;

  REQUIRE(ms1 == ms2);

  Vector d1(6,12);
  Vector d3(6,3);

  REQUIRE(d1 / 4 == d3);

  d1 /= 4;

  REQUIRE(d1 == d3);

  Vector v1(5,4);
  Matrix mt1(5,3.0);
  Vector va2(5,60);

  REQUIRE(v1 * mt1 == va2);

  v1 *= mt1;

  REQUIRE(v1 == va2);


}

TEST_CASE("Matrix"){

  Matrix mat2(4,7.0);
  Matrix mat3(4,2.5);
   Matrix mat4(4,4.5);

   REQUIRE(mat2 - mat3 == mat4);
  mat2 -= mat3;
  REQUIRE(mat2 == mat4);
  REQUIRE(mat3 != mat4);

  REQUIRE((mat3 + mat2).get(0,0) == 7.0);

  mat3 += mat2;

  REQUIRE(mat3.get(2,2) == 7.0);


  Matrix matr(3,4,2);
  matr.reshape(2,6);
  REQUIRE(matr.shape().first == 2);
  REQUIRE(matr.shape().second == 6);

 Matrix matr2(3,7.0);
  Matrix matr3(3,8.0);
  Matrix matr4(3,168.0);

  REQUIRE(matr2 * matr3 == matr4);

  matr2 *= matr3;

 REQUIRE(matr2 == matr4);

  Matrix matr5(6,2.0);
  Matrix matr6(6,14.0);

  Matrix matr8(3,7.0);
  Vector vctr(3,8);
  Vector vctr1(3,168);

  REQUIRE(matr8 * vctr == vctr1);

  Matrix matr10(8,1,3);
  Matrix matr11(1,8,3);
  matr10.transpose();

  REQUIRE(matr10 == matr11);

  Matrix matr20(1,5.0);
  REQUIRE(matr20.det() == 5);

}
