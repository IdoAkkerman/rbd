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
#ifndef RBD_VECTOR_HPP
#define RBD_VECTOR_HPP

#include <iostream>
#include <json/json.h>
#include <assert.h>

//=========================================================
// Define basic Vector object
//=========================================================
typedef std::vector<double> Vector;

//=========================================================
// Vector resize
//=========================================================
void Resize(int dim, Vector&v)
{
   v.resize(dim);
   for(double &d: v)
      d = 0;
}

//=========================================================
// Vector-Scalar product
//=========================================================
Vector operator*(const double &c,const Vector &v)
{
   int n = v.size();

   Vector prod;
   Resize(n,prod);

   for(int i = 0; i < n; i++)
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
// Vector-Vector innerprod
//=========================================================
double operator*(const Vector &v1, const Vector &v2)
{
   assert(v1.size() == v2.size());
   int n = v1.size();

   double prod = 0.0;
   for(int i = 0; i < n; i++)
   {
      prod += v1[i] * v2[i];
   }
   return prod;
}

//=========================================================
// Vector-Vector sum
//=========================================================
Vector operator+(const Vector &v1, const Vector &v2)
{
   assert(v1.size() == v2.size());
   int n = v1.size();

   Vector sum;
   Resize(n,sum);

   for(int i = 0; i < n; i++)
   {
      sum[i] = v1[i] + v2[i];
   }
   return sum;
}

//=========================================================
// Vector-Vector dif
//=========================================================
Vector operator-(const Vector &v1, const Vector &v2)
{
   assert(v1.size() == v2.size());
   int n = v1.size();

   Vector dif;
   Resize(n,dif);

   for(int i = 0; i < n; i++)
   {
      dif[i] = v1[i] - v2[i];
   }
   return dif;
}

//=========================================================
// Squared norm
//=========================================================
double Norm(const Vector &v)
{
   double norm = 0.0;
   for (double d: v)
      norm += d*d;
   return norm;
}

//=========================================================
// Vector print routine
//=========================================================
void PrintVector(std::ostream &out, Vector &vec)
{
   out<<"[";
   for(double d: vec)
      out<<d<<" ";
   out<<"]\n";
}

//=========================================================
// Vector Json conversion routines
//=========================================================
void json2vector(Json::Value& jv, Vector &v)
{
   int dim = jv.size();
   v.resize(dim);
   for (int i = 0; i < dim; ++i)
      v[i] = jv[i].asDouble();
}

void json2vector(Json::Value& jv, const char *key, Vector &v)
{
   int dim = v.size();
   if (!jv.isMember(key)) return;
   json2vector(jv["x"], v);
   assert(dim == v.size());
}

#endif
