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
// Matrix-Matrix-Matrix product for tensor rotation
//=========================================================
void RIRT(const Matrix &R, const Matrix &I0, Matrix &I)
{
   assert(R.size() == I0.size());  // Assume matrices are square
   assert(R.size() == I.size());   // Assume matrices are square
   int n = R.size();
   int i,j,k,l;
   double tmp;
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         I[i][j] = 0.0;
         for(k = 0; k< n; k++)
         {
            tmp = 0.0;
            for(l = 0; l < n; l++)
               tmp += R[j][l]*I0[k][l];
            I[i][j] += R[i][k]*tmp;
         }
      }
   }
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
