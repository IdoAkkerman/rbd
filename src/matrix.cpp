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

#include <iostream>
#include <json/json.h>
#include <assert.h>
#include <cmath>

#include "matrix.hpp"

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Matrix&A)
{
   A.resize(dim);
   for (Vector &vec: A)
   {
      Resize(dim, vec);
   }
}

//=========================================================
// Check Matrix properties
//=========================================================
bool isSymmetric(const Matrix&A)
{
   int n = A.size();
   bool sym = true;
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
         if (fabs(A[j][i] - A[i][j]) > 1e-12) { sym = false; }
   }
   return sym;
}

bool isDiagionallyDominant(const Matrix&A)
{

   int n = A.size();
   bool pd = true;
   for (int i = 0; i < n; i++)
   {
      double sum = 0;
      for (int j = 0; j < n; j++)
         if (j != i) { sum += fabs(A[i][j]); }

      std::cout<<    fabs(A[i][i])<<" < "<<sum<<std::endl;
      if (fabs(A[i][i]) < sum) { pd = false; }
   }
   return pd;
}

double det(const Matrix &A)
{
   assert((A.size() == 2) || (A.size() == 3));
   double det;
   if (A.size() == 2)
   {
      det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
   }
   else if (A.size() == 3)
   {
      det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
            A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
            A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
   }
   else
   {
      abort();
   }
   return det;
}

bool isInvertible(const Matrix&A)
{
   return (fabs(det(A)) > 1e-12);
}

double trace(const Matrix &A)
{
   int n = A.size();
   double tr = 0.0;
   for (int i = 0; i < n; i++)
   {
      tr += A[i][i];
   }

   return tr;
}

bool isPositiveDefinite(const Matrix&A)
{
   assert((A.size() == 2) || (A.size() == 3));

   if (A.size() == 2)
   {
      return (trace(A) > 1e-12) && (det(A) > 1e-12);
   }
   else if (A.size() == 3)
   {
      double I2 = A[0][0] * A[1][1]
                  + A[1][1] * A[2][2]
                  + A[0][0] * A[2][2]
                  - A[0][1] * A[1][0]
                  - A[1][2] * A[2][1]
                  - A[0][2] * A[2][0];
      return (I2 > 1e-12) && (trace(A) > 1e-12) &&(det(A) > 1e-12);
   }
   else
   {
      abort();
   }
   return false;
}

bool isSPD(const Matrix&A)
{
   return isSymmetric(A) && isPositiveDefinite(A);
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
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         At[j][i]= A[i][j];
      }
   }
   return At;
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

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
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
// Matrix-Vector product
//=========================================================
Vector operator*(const Matrix &A, const Vector &x)
{
   int i,j;
   int m = A.size();
   int n = x.size();

   Vector prod(m);

   for (i = 0; i < m; i++)
   {
      prod[i] = 0.;
      for (j = 0; j < n; j++)
      {
         prod[i] += A[i][j]*x[j];
      }
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

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         prod[i][j] = 0.;
         for (k = 0; k< n; k++)
         {
            prod[i][j] += A[i][k]*B[k][j];
         }
      }
   }
   return prod;
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

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         sum[i][j] = A[i][j] + B[i][j];
      }
   }
   return sum;
}


//=========================================================
// Matrix-Matrix dif
//=========================================================
Matrix operator-(const Matrix &A, const Matrix &B)
{
   assert(A.size() == B.size());  // Assume matrices are square
   int n = A.size();

   int i,j,k;

   Matrix dif;
   Resize(n,dif);

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         dif[i][j] = A[i][j] - B[i][j];
      }
   }
   return dif;
}


//=========================================================
// Squared norms
//=========================================================
double Norm(const Matrix &A)
{
   double norm = 0.0;
   for (Vector vec: A)
   {
      norm += Norm(vec);
   }
   return norm;
}

//=========================================================
// Compute Inverse of Matrix
//=========================================================
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

//=========================================================
// Generate the skew Matrix from a vector
//=========================================================
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
// Matrix print routine
//=========================================================
void PrintMatrix(std::ostream &out, Matrix &mat)
{
   out<<"\n";
   for (Vector &vec: mat)
   {
      PrintVector(out, vec);
   }
}

//=========================================================
// Matrix Json conversion routines
//=========================================================
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

