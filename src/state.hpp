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

//=========================================================
// Define the rigid body class
//=========================================================
class State
{
public:
   int dim;
   Vector x,v;  // Displacement & Velocity;
   Matrix R;    // Rotation matrix
   Vector w;    // Angular velocity

   Matrix I;    // Inertia tensor --> in current configuration

   State() {};
   State(State &org);

   void Initialize(int dim_);
   void SetInertiaTensor(Matrix &I0);
   void CheckRotation();
   void Read(Json::Value &data);
   void Print(std::ostream &out);
};

#endif
