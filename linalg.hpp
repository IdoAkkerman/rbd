/*
* RBD Rigid body dynamics adaptor for precice
* Copyright (C) 2025  I.Akkerman@tudelft.nl
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef RBD_LINALG_HPP
#define RBD_LINALG_HPP

#include <iostream>
#include <json/json.h>
#include <assert.h>

//=========================================================
// Define basic Vector and Matrix objects
//=========================================================
typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Vector&v)
{
   v.resize(dim);
   for(double &d: v)
      d = 0;
}

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Matrix&A)
{
   A.resize(dim);
   for(Vector &vec: A)
      Resize(dim, vec);
}

//=========================================================
// Matrix transpose
//=========================================================
Matrix Transpose(Matrix&A)
{
   int i,j;
   int n = A.size();
   Matrix At;
   Resize(n,At);
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
          At[j][i]= A[i][j];
   }
   return At;
}

//=========================================================
// Matrix-Vector product
//=========================================================
Vector operator*(const Matrix &A, const Vector &x)
{
   int i,j;
   int m = A.size();
   int n = x.size();

   Vector prod(m);

   for(i = 0; i < m; i++)
   {
      prod[i] = 0.;
      for(j = 0; j < n; j++)
          prod[i] += A[i][j]*x[j];
   }
   return prod;
}

//=========================================================
// Matrix-Matrix product
//=========================================================
Matrix operator*(const Matrix &A, const Matrix &B)
{
   assert(A.size() == B.size());  // Assume matrices are square
   int n = A.size();

   int i,j,k;

   Matrix prod;
   Resize(n,prod);

   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         prod[i][j] = 0.;
         for(k = 0; k< n; k++)
            prod[i][j] += A[i][k]*B[k][j];
      }
   }
   return prod;
}

//=========================================================
// Matrix-Scalar product
//=========================================================
Matrix operator*(const double &c,const Matrix &A)
{
   // Assume matrices are square
   int n = A.size();

   int i,j,k;

   Matrix prod;
   Resize(n,prod);

   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
            prod[i][j] = A[i][j]*c;
      }
   }
   return prod;
}

Matrix operator*(const Matrix &A, const double &c)
{
   return c*A;
}

//=========================================================
// Vector-Scalar product
//=========================================================
Vector operator*(const double &c,const Vector &v)
{
   int n = v.size();
   int i;

   Vector prod;
   Resize(n,prod);

   for(i = 0; i < n; i++)
   {
      prod[i] = v[i]*c;
   }
   return prod;
}

Vector operator*(const Vector &v,const double &c)
{
   return c*v;
}

//=========================================================
// Matrix-Matrix sum
//=========================================================
Matrix operator+(const Matrix &A, const Matrix &B)
{
   assert(A.size() == B.size());  // Assume matrices are square
   int n = A.size();

   int i,j,k;

   Matrix sum;
   Resize(n,sum);

   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         sum[i][j] = A[i][j] + B[i][j];
      }
   }
   return sum;
}

//=========================================================
// Vector-Vector sum
//=========================================================
Vector operator+(const Vector &v1, const Vector &v2)
{
   assert(v1.size() == v2.size());
   int n = v1.size();

   int i;

   Vector sum;
   Resize(n,sum);

   for(i = 0; i < n; i++)
   {
      sum[i] = v1[i] + v2[i];
   }
   return sum;
}


//=========================================================
// Squared norms
//=========================================================
double Norm(const Vector &v)
{
   double norm = 0.0;
   for (double d: v)
      norm += d*d;
   return norm;
}

double Norm(const Matrix &A)
{
   double norm = 0.0;
   for (Vector vec: A)
      norm += Norm(vec);
   return norm;
}


Matrix Inverse(const Matrix &m)
{
   assert(m.size() == 3);
   double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

   double invdet = 1 / det;

   Matrix inv;
   Resize(3,  inv);
   inv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
   inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
   inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
   inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
   inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
   inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
   inv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
   inv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
   inv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

   return inv;
}



Matrix Skew(const Vector &v)
{
   int n =v.size();
   Matrix sk;
   Resize(n, sk);
   assert(n == 3);  // assumes 3D
   sk[0][1] = -v[2];
   sk[0][2] = v[1];
   sk[1][2] = -v[0];

   sk[1][0] = v[2];
   sk[2][0] = -v[1];
   sk[2][1] = v[0];
   return sk;
}

//=========================================================
// Matrix and Vector print routines
//=========================================================
void Print(std::ostream &out, Vector &vec)
{
   out<<"[";
   for(double d: vec)
      out<<d<<" ";
   out<<"]\n";
}

void Print(std::ostream &out, Matrix &mat)
{
   out<<"\n";
   for(Vector &vec: mat)
      Print(out, vec);
}

//=========================================================
// Matrix and Vector Json conversion routines
//=========================================================

void json2vector(Json::Value& jv, Vector &v)
{
   int dim = jv.size();
   v.resize(dim);
   for (int i = 0; i < dim; ++i)
      v[i] = jv[i].asDouble();
}

void json2vector(Json::Value& jv, char *key, Vector &v)
{
   int dim = v.size();
   if (!jv.isMember("key")) return;
   json2vector(jv["x"], v);
   assert(dim == v.size());
}

void json2matrix(Json::Value& jv, Matrix &m)
{
   int dim = jv.size();
   m.resize(dim);
   for (int i = 0; i < dim; ++i)
   {
      json2vector(jv[i],m[i]);
      assert(m[i].size() == dim);   // Force the matrix to be square
   }
}

#endif
