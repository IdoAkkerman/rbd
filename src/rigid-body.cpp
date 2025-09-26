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

#include <assert.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
//#include <json/json.h>

#include "vector.hpp"
#include "matrix.hpp"
#include "state.hpp"
#include "rigid-body.hpp"

// Defacto Constructor
void RigidBody::Construct(Json::Value &data)
{
   // Read object data
   dim = data["dimension"].asInt();
   m = data["mass"].asDouble();
   json2matrix(data["I0"],I0);
   PrintMatrix(std::cout, I0);
   assert(I0.size() == dim);
   assert(isSPD(I0));

   // Read integration data
   dt_max = data.get("dt", 9999999.0).asDouble();

   // Initialize states
   state_new.Initialize(dim);
   state_old.Initialize(dim);

   if (data.isMember("initial_condition"))
   {
      state_old.Read(data["initial_condition"]);
   }
   state_old.SetInertiaTensor(I0);
}

// Constructor
RigidBody::RigidBody(std::string fileName, std::string rbName)
{
   Json::Value data;
   std::ifstream file(fileName.c_str(), std::ifstream::binary);
   if (!file)
   {
      std::cout<<"File "<<fileName<<" could not be opened!\n";
      abort();
   }
   file >> data;
   Construct(data[rbName]);
}

// Constructor
RigidBody::RigidBody(std::string fileName)
{
   std::ifstream file(fileName.c_str(), std::ifstream::binary);
   if (!file)
   {
      std::cout<<"File "<<fileName<<" could not be opened!\n";
      abort();
   }
   Json::Value data;
   file >> data;
   Construct(data);
}

// Constructor
RigidBody::RigidBody(Json::Value &data)
{
   Construct(data);
}

// Print
void RigidBody::Print(std::ostream &out)
{
   // Object data
   out<<"Dimension = "<<dim<<std::endl;
   out<<"Mass = "<<m<<std::endl;
   out<<"I0 = "; PrintMatrix(out, I0);

   // State data
   state_old.Print(out);
   state_new.Print(out);
}

// Time integration routine
void RigidBody::computeMotion(double dt, Vector &forces,
                              int iter_max, double tol)
{
   Vector F(dim);
   Vector M(dim);
   int i = 0;
   for (double &d: F)
   {
      d = forces[i++];
   }
   for (double &d: M)
   {
      d = forces[i++];
   }

   state_new.v =  state_old.v + (dt/m)*F;
   state_new.x =  state_old.x + (dt/2)*(state_new.v+state_old.v);

   Matrix eye;
   Resize(dim, eye);
   for (int i = 0; Vector & vec: eye)
   {
      vec[i++] = 1.0;
   }

   int it;
   for (it = 0; it < iter_max; it++)
   {
      state_new.SetInertiaTensor(I0);
      Matrix I_inv = Inverse(state_new.I);
      state_new.w =  I_inv*(state_old.I*state_old.w + dt*I_inv*M);

      Matrix Q = Skew(0.5*(state_new.w+state_old.w));
      Matrix Qp = eye + (dt/2)*Q;
      Matrix Qm_inv = Inverse(eye + (-dt/2)*Q);
      Matrix Rn= Qm_inv*Qp*state_old.R;
      double norm = Norm(state_new.R - Rn);
      state_new.R =  Rn;

      if (norm < tol) { break; }
   }

   std::cout<<"iterations = "<<it<<std::endl;
}
