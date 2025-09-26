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
#ifndef RBD_MATRIX_HPP
#define RBD_MATRIX_HPP

#include <iostream>
#include <json/json.h>
#include <assert.h>
#include <cmath>

#include "vector.hpp"

//=========================================================
// Define basic Matrix object
//=========================================================
typedef std::vector<Vector> Matrix;

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Matrix&A);

//=========================================================
// Check Matrix properties
//=========================================================
double det(const Matrix &A);
double trace(const Matrix &A);

bool isSymmetric(const Matrix&A, double tol = 1e-12);
bool isDiagionallyDominant(const Matrix&A);
bool isInvertible(const Matrix&A, double tol = 1e-12);
bool isPositiveDefinite(const Matrix&A, double tol = 1e-12);
bool isSPD(const Matrix&A, double tol = 1e-12);

//=========================================================
// Matrix transpose
//=========================================================
Matrix Transpose(Matrix&A);

//=========================================================
// Matrix-Scalar product
//=========================================================
Matrix operator*(const double &c,const Matrix &A);
Matrix operator*(const Matrix &A, const double &c);

//=========================================================
// Matrix-Vector product
//=========================================================
Vector operator*(const Matrix &A, const Vector &x);

//=========================================================
// Matrix-Matrix product
//=========================================================
Matrix operator*(const Matrix &A, const Matrix &B);

//=========================================================
// Matrix-Matrix sum
//=========================================================
Matrix operator+(const Matrix &A, const Matrix &B);

//=========================================================
// Matrix-Matrix dif
//=========================================================
Matrix operator-(const Matrix &A, const Matrix &B);

//=========================================================
// Squared norms
//=========================================================
double Norm(const Matrix &A);

//=========================================================
// Compute Inverse of Matrix
//=========================================================
Matrix Inverse(const Matrix &m);

//=========================================================
// Generate the skew Matrix from a vector
//=========================================================
Matrix Skew(const Vector &v);

//=========================================================
// Matrix print routine
//=========================================================
void PrintMatrix(std::ostream &out, Matrix &mat);

//=========================================================
// Matrix Json conversion routines
//=========================================================
void json2matrix(Json::Value& jv, Matrix &m);

#endif
