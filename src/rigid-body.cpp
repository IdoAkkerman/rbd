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
   rdim = (dim*(dim-1))/2;

   // Translational parameters
   m = data["mass"].asDouble();

   Resize(dim, Ct);
   json2matrix(data,"Ct",Ct);

   Resize(dim, Kt);
   json2matrix(data,"Kt",Kt);

   // Rotational parameters
   json2matrix(data["I0"],I0);
   PrintMatrix(std::cout, I0);
   assert(I0.size() == rdim);
   assert(isSPD(I0));

   Resize(rdim, Cr);
   json2matrix(data,"Cr",Cr);

   Resize(dim, Kr);
   json2matrix(data,"Kr",Kr);

   // Integration parameters
   dt_max = data.get("dt", 9999999.0).asDouble();
   theta_max = data.get("theta_max", 0.5).asDouble();

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
   std::ifstream file(fileName.c_str(), std::ifstream::binary);
   if (!file)
   {
      std::cout<<"File "<<fileName<<" could not be opened!\n";
      abort();
   }
   Json::Value data;
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
   PrintMatrix(out, Ct, "Ct");
   PrintMatrix(out, Kt, "Kt");

   PrintMatrix(out, I0, "I0");
   PrintMatrix(out, Cr, "Cr");
   PrintMatrix(out, Kr, "Kr");
}

double RigidBody::GetTimeStep()
{
   return std::fmin(dt_max, theta_max/sqrt(Norm(state_old.w) + 1e-24));
};

// Time integration routine
void RigidBody::ComputeMotion(double dt, Vector &F, Vector &M,
                              int iter_max, double tol)
{
   state_new.t =  state_old.t + dt;

   Vector F2(2*dim);

   for (int i = 0; i < dim; i++)
   {
      F2[i] = F[i];
   }

   Matrix eye;
   Resize(dim, eye);
   for (int i = 0; Vector &vec: eye)
   {
      vec[i++] = 1;
   }

   if (dt != dt_A)
   {
      dt_A = dt;

      Matrix Mt;
      Resize(dim, Mt);
      for (int i = 0; Vector &vec: Mt)
      {
         vec[i++] = m;
      }

      Matrix A1;
      Resize(2*dim, A1);
      SetBlock(A1,0,0,Mt + 0.5*dt*Ct);
      SetBlock(A1,0,1, 0.5*dt*Kt);
      SetBlock(A1,1,0,-0.5*dt*eye);
      SetBlock(A1,1,1,eye);

      Resize(2*dim, A0);
      SetBlock(A0,0,0,Mt - 0.5*dt*Ct);
      SetBlock(A0,0,1,-0.5*dt*Kt);
      SetBlock(A0,1,0, 0.5*dt*eye);
      SetBlock(A0,1,1,eye);

      A1_inv = Inverse(A1);
   }

   Vector xv0(2*dim);
   for (int j = 0; j < dim; j++)
   {
      xv0[j] = state_old.v[j];
      xv0[j + dim] = state_old.x[j];
   }

   Vector xv = A1_inv*(F2 + A0*xv0);
   for (int j = 0; j < dim; j++)
   {
      state_new.v[j] = xv[j];
      state_new.x[j] = xv[j + dim];
   }

   for (iter = 0; iter < iter_max; iter++)
   {
      state_new.SetInertiaTensor(I0);
      Matrix I_inv = Inverse(state_new.I + (0.5*dt)*Cr);
      Vector M_mod = M;
      if (hasKr)
      {
         state_new.theta = GetRotation(state_new.R);
         M_mod = M - 0.5*Kr*(state_old.theta + state_new.theta);
      }
      Vector test = state_old.w + dt*M_mod;

      state_new.w =  I_inv*((state_old.I - (0.5*dt)*Cr)*state_old.w + dt*M_mod);

      Matrix Q = Skew(0.5*(state_new.w+state_old.w));
      Matrix Qp = eye + (dt/2)*Q;
      Matrix Qm_inv = Inverse(eye - (dt/2)*Q);

      Matrix Rn= Qm_inv*Qp*state_old.R;

      double norm = Norm(state_new.R - Rn);
      state_new.R =  Rn;

      if (norm < tol) { break; }
   }
   if (hasKr) { state_new.theta = GetRotation(state_new.R); }
   if (verbose)
   {
      std::cout<<"iter = "<<iter<<std::endl;
      std::cout<<"dtheta = "<<dt*sqrt(Norm(state_new.w)) <<std::endl;
   }
}
