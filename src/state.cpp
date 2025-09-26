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
#include <cmath>

#include "state.hpp"

//=========================================================
// Convert the angles to a rotation matrix
// Using the Tait-Bryan version of the Euler-angels
// Specifically the nautical or Cardan angles
//=========================================================
Matrix GetRotation(Vector rad)
{
   Matrix R;
   if (rad.size() == 1)
   {
      Resize(2,R);
      R[0][0] = cos(rad[0]);
      R[0][1] = -sin(rad[0]);
      R[1][0] = -R[0][1];
      R[1][1] = R[0][0];
   }
   else if (rad.size() == 3)
   {
      Matrix Zalpha, Ybeta, Xgamma;
      Resize(3, Zalpha);
      Zalpha[0][0] = cos(rad[0]);
      Zalpha[0][1] = -sin(rad[0]);
      Zalpha[1][0] = -Zalpha[0][1];
      Zalpha[1][1] = Zalpha[0][0];
      Zalpha[2][2] = 1.0;

      Resize(3, Ybeta );
      Ybeta[0][0] = cos(rad[1]);
      Ybeta[0][2] = sin(rad[1]);  // Flipped!!
      Ybeta[1][1] = 1.0;
      Ybeta[2][0] = -Ybeta[0][2];
      Ybeta[2][2] = Ybeta[0][0];

      Resize(3, Xgamma);
      Xgamma[0][0] = 1.0;
      Xgamma[1][1] = cos(rad[2]);
      Xgamma[1][2] = -sin(rad[2]);
      Xgamma[2][1] = -Xgamma[1][2];
      Xgamma[2][2] = Xgamma[1][1];

      R = Zalpha*Ybeta*Xgamma;
   }
   else
   {
      std::cout<<"Wrong number of angles specified";
      abort();
   }
   return R;
}

//=========================================================
// Define the rigid body class
//=========================================================
State::State(State &org)
{
   dim = org.dim;
   x = org.x;
   v = org.v;
   R = org.R;
   w = org.w;
   I = org.I;
}

void State::Initialize(int dim_)
{
   dim = dim_;
   Resize(dim, x);
   Resize(dim, v);
   Resize((dim*(dim-1))/2, w);
   Resize(dim, R);
   for (int i = 0; Vector &vec: R)
   {
      vec[i++] = 1.0;
   }
   Resize(dim, I);
}

void State::SetInertiaTensor(Matrix &I0)
{
   assert(I0.size() == (dim*(dim-1))/2);

   if (dim == 2)
   {
      I = I0;
   }
   else if (dim == 3)
   {
      I = R*I0*Transpose(R);
   }
   else
   {
      abort();
   }
}

void State::CheckRotation(double tol)
{
   // Check rotation matrix
   Matrix id = R*Transpose(R);

   for (int i = 0; Vector &vec: id)
   {
      vec[i++] -= 1.0;
   }
   double norm = Norm(id);
   if (norm > tol)
   {
      std::cout<<" Norm = " << norm<<std::endl;
      std::cout<<" R = "; PrintMatrix(std::cout, R);
      Matrix id = R*Transpose(R);
      std::cout<<" R*R^t = "; PrintMatrix(std::cout, id);
   }
   assert(norm < tol);
}

void State::Read(Json::Value &data)
{
   // Read state data -- zero if not defined
   // Displacement
   json2vector(data, "x", x);

   // Velocity
   json2vector(data, "v", v);

   // Angular Velocity
   json2vector(data, "w", w);

   // Rotation
   if (data.isMember("angles_deg"))
   {
      Vector angles_deg;
      json2vector(data["angles_deg"], angles_deg);
      assert(angles_deg.size() == (dim*(dim-1))/2);
      for (double &d: angles_deg)
      {
         d *= M_PI/180;
      }

      R = GetRotation(angles_deg);
   }
   else if (data.isMember("radians"))
   {
      Vector radians;
      json2vector(data["radians"], radians);
      assert(radians.size()  == (dim*(dim-1))/2);
      R = GetRotation(radians);
   }
   else if (data.isMember("R"))
   {
      json2matrix(data["R"], R);
      assert(R.size() == dim);
   }
   else
   {
      Resize(dim, R);
      for (int i = 0; Vector &vec: R)
      {
         vec[i++] = 1.0;
      }
   }
};

void State::Print(std::ostream &out)
{
   out<<"x = "; PrintVector(out, x);
   out<<"v = "; PrintVector(out, v);

   out<<"R = "; PrintMatrix(out, R);
   out<<"w = "; PrintVector(out, w);
}


