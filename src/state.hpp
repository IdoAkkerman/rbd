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
#ifndef RBD_STATE_HPP
#define RBD_STATE_HPP

#include <json/json.h>
#include <iostream>

#include "vector.hpp"
#include "matrix.hpp"


Matrix GetRotation(const Vector &rad);
Vector GetRotation(const Matrix &R);

//=========================================================
// Define the rigid body class
//=========================================================
class State
{
public:
   int dim, rdim;
   double t = -1.0;
   Vector x,v;     // Displacement & Velocity;
   Matrix R;       // Rotation matrix
   Vector theta;   // Rotation angels
   Vector w;       // Angular velocity
   Matrix I;       // Inertia tensor --> in current configuration

   State() {};
   State(int dim);
   State(const State &org);

   void Initialize(int dim_);
   void SetInertiaTensor(Matrix &I0);
   void CheckRotation(double tol = 1e-12) const;
   bool Equal(const State &org, double tol = 1e-12) const;
   bool Read(std::string fileName, std::string rbName);
   void Read(Json::Value &data);
   void PrettyPrint(std::ostream &out = std::cout) const;
   void PrintJson(std::ostream &out = std::cout) const;
   void PrintHeader(std::ostream &out) const;
   void Print(double dt, std::ostream &out) const;
};

#endif
