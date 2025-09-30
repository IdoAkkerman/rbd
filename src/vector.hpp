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
void Resize(int dim, Vector&v);

//=========================================================
// Vector-Scalar product
//=========================================================
Vector operator*(const double &c,const Vector &v);
Vector operator*(const Vector &v,const double &c);

//=========================================================
// Vector-Vector innerprod
//=========================================================
double operator*(const Vector &v1, const Vector &v2);

//=========================================================
// Vector-Vector sum
//=========================================================
Vector operator+(const Vector &v1, const Vector &v2);

//=========================================================
// Vector-Vector dif
//=========================================================
Vector operator-(const Vector &v1, const Vector &v2);

//=========================================================
// Squared norm
//=========================================================
double Norm(const Vector &v);

//=========================================================
// Vector print routine
//=========================================================
void PrintVector(std::ostream &out, const Vector &vec,
                 const std::string &key = "");

//=========================================================
// Vector Json conversion routines
//=========================================================
void json2vector(Json::Value& jv, Vector &v);
void json2vector(Json::Value& jv, const char *key, Vector &v);

Json::Value to_json(const Vector &v);

#endif
