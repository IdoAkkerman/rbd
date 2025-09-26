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
#ifndef RBD_RIGID_BODY_HPP
#define RBD_RIGID_BODY_HPP

#include <assert.h>
#include <iostream>
#include <cmath>
#include <json/json.h>

#include "vector.hpp"
#include "matrix.hpp"
#include "state.hpp"

//=========================================================
// Define the rigid body class
//=========================================================
class RigidBody
{
private:
   int dim = -1; // Dimension of the problem (2D/3D)

   double m;     // Mass
   Matrix I0;    // Mass moment of inertia tensor

   double dt_max;
   int iter = -1;

   State state_old, state_new;

   // Defacto contructor
   void Construct(Json::Value &data);

public:
   // Contructors
   RigidBody(std::string fileName, std::string rbName);
   RigidBody(std::string fileName);
   RigidBody(Json::Value &data);

   // Info
   int GetDimension() {return dim;};
   int GetIterationCount() {return iter;};
   void Print(std::ostream &out = std::cout);
   State& GetNewState() { return state_new; };
   State& GetOldState() { return state_old; };

   // Time stepping routines
   double GetTimeStep() {return dt_max;};
   void ComputeMotion(double dt, Vector &forces,
                      int iter_max = 100,
                      double tol = 1e-8);

   void ReloadOldState() {}; // Not needed
   void SaveOldState() {};  // Not needed
   void EndTimeStep() { state_old = state_new; };
};

#endif
