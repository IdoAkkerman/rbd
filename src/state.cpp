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
#include <fstream>
#include <json/json.h>

#include "state.hpp"

//=========================================================
// Convert the angles to a rotation matrix
// Using the Tait-Bryan version of the Euler-angels
// Specifically the nautical or Cardan angles
//=========================================================
Matrix GetRotation(const Vector &rad)
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
      Matrix Zgamma, Ybeta, Xalpha;

      Resize(3, Xalpha);
      Xalpha[0][0] = 1.0;
      Xalpha[1][1] = cos(rad[0]);
      Xalpha[1][2] = -sin(rad[0]);
      Xalpha[2][1] = -Xalpha[1][2];
      Xalpha[2][2] = Xalpha[1][1];

      Resize(3, Ybeta );
      Ybeta[0][0] = cos(rad[1]);
      Ybeta[0][2] = sin(rad[1]);  // Flipped!!
      Ybeta[1][1] = 1.0;
      Ybeta[2][0] = -Ybeta[0][2];
      Ybeta[2][2] = Ybeta[0][0];

      Resize(3, Zgamma);
      Zgamma[0][0] = cos(rad[2]);
      Zgamma[0][1] = -sin(rad[2]);
      Zgamma[1][0] = -Zgamma[0][1];
      Zgamma[1][1] = Zgamma[0][0];
      Zgamma[2][2] = 1.0;

      R = Xalpha*Ybeta*Zgamma;
   }
   else
   {
      std::cout<<"Wrong number of angles specified";
      abort();
   }
   return R;
}

Vector GetRotation(const Matrix &R)
{
   int dim = R.size();
   Vector theta((dim*(dim-1))/2);
   if (R.size() == 2)
   {
      theta[0] = atan(-R[1][0]/R[1][1]);
   }
   else if (R.size()== 3)
   {
      theta[0] = atan2(-R[1][2],R[2][2]);
      theta[1] = asin(R[0][2]);
      theta[2] = atan2(-R[0][1],R[0][0]);
      if (fabs(theta[1]) > 1.5) { std::cout<<" "<<std::endl; }
   }
   else
   {
      abort();
   }
   return theta;
}


//=========================================================
// Define the rigid body class
//=========================================================
State::State(int dim_)
{
   Initialize(dim_);
}

State::State(const State &org)
{
   dim = org.dim;
   rdim = org.rdim;
   t = org.t;
   x = org.x;
   v = org.v;
   theta = org.theta;
   R = org.R;
   w = org.w;
   I = org.I;
}

void State::Initialize(int dim_)
{
   dim = dim_;
   rdim = (dim*(dim-1))/2;
   t = 0.0;
   Resize(dim, x);
   Resize(dim, v);
   Resize(rdim, w);
   Resize(rdim, theta);
   Resize(dim, R);
   for (int i = 0; Vector &vec: R)
   {
      vec[i++] = 1.0;
   }
   Resize(dim, I);
}

void State::SetInertiaTensor(Matrix &I0)
{
   assert(I0.size() == rdim);

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

void State::CheckRotation(double tol) const
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

bool State::Equal(const State &org, double tol) const
{
   assert(dim == org.dim);
   assert(rdim == org.rdim);
   double dt = t - org.t;
   Vector dx = x - org.x;
   Vector dv = v - org.v;
   Vector dtheta = theta - org.theta;
   Matrix dR = R - org.R;
   Vector dw = w - org.w;

   bool check = (fabs(dt) < tol) &&
                (Norm(dx) < tol) &&
                (Norm(dv) < tol) &&
                (Norm(dtheta) < tol) &&
                (Norm(dR) < tol) &&
                (Norm(dw) < tol);

   if (!check)
   {
      std::cout<<"dt = "<<dt<<std::endl;
      PrintVector(std::cout, dx, "dx");
      PrintVector(std::cout, dv, "dv");
      PrintVector(std::cout, dtheta, "dtheta");
      PrintMatrix(std::cout, dR, "dR");
      PrintVector(std::cout, dw, "dw");
   }

   return check;
}

bool State::Read(std::string fileName, std::string stateName)
{
   std::ifstream file(fileName.c_str(), std::ifstream::binary);
   if (!file)
   {
      std::cout<<"File "<<fileName<<" could not be opened!\n";
      return false;
   }
   Json::Value data;
   file >> data;
   if (!data.isMember(stateName))
   {
      std::cout<<"Could not find "<<stateName<<" in file "<<fileName<<"\n";
      return false;
   }
   Read(data[stateName]);
   return true;
}

void State::Read(Json::Value &data)
{
   // Read state data -- zero if not defined
   t = data["t"].asDouble();

   // Displacement
   json2vector(data, "x", x);

   // Velocity
   json2vector(data, "v", v);

   // Angular Velocity
   json2vector(data, "w", w);

   // Rotation matrix
   if (data.isMember("angles_deg"))
   {
      json2vector(data["angles_deg"], theta);
      assert(theta.size() == rdim);
      for (double &d: theta)
      {
         d *= M_PI/180;
      }

      R = GetRotation(theta);
   }
   else if (data.isMember("radians"))
   {
      json2vector(data["radians"], theta);
      assert(theta.size()  == rdim);
      R = GetRotation(theta);
   }
   else if (data.isMember("R"))
   {
      json2matrix(data["R"], R);
      assert(R.size() == dim);
      theta = GetRotation(R);
   }
   else
   {
      Resize(dim, theta);
      Resize(dim, R);
      for (int i = 0; Vector &vec: R)
      {
         vec[i++] = 1.0;
      }
   }

   // Rotation angles
   json2vector(data, "theta", theta);
};

void State::PrettyPrint(std::ostream &out) const
{
   out<<"t = "<<t<<std::endl;
   PrintVector(out, x, "x");
   PrintVector(out, v, "v");

   PrintVector(out, theta, "theta");
   PrintMatrix(out, R, "R");
   PrintVector(out, w, "w");
}

void State::PrintJson(std::ostream &out) const
{
   Json::Value j(Json::objectValue);
   j["t"] = t;

   //  Json::json j;// = x;

   j["x"] = to_json(x);
   j["v"] = to_json(v);

   j["theta"] = to_json(theta);
   j["R"] = to_json(R);
   j["w"] = to_json(w);

   out<<j<<std::endl;
}

void State::PrintHeader( std::ostream &out) const
{
   int ii = 0;
   out<<"#"<<ii++<<":t "<<ii++<<":dt ";
   for (int i = 0; const double &d: x) { out<<ii++<<":x_"<<i++<<" "; }
   for (int i = 0; const double &d: v) { out<<ii++<<":v_"<<i++<<" "; }

   for (int i = 0; const double &d: theta) { out<<ii++<<":theta_"<<i++<<" "; }
   for (int i = 0; const Vector &vec: R)
   {
      for (int j = 0; const double &d: vec) { out<<ii++<<":R"<<i<<j++<<" "; }
   }
   for (int i = 0; const double &d: w) { out<<ii++<<":w_"<<i++<<" "; }
   out<<std::endl;
}

void State::Print(double dt, std::ostream &out) const
{
   out<<t<<" "<<dt;
   for (const double &d: x) { out<<" "<<d; }
   for (const double &d: v) { out<<" "<<d; }

   for (const double &d: theta) { out<<" "<<d; }
   for (const Vector &vec: R)
      for (const double &d: vec) { out<<" "<<d; }
   for (const double &d: w) { out<<" "<<d; }
   out<<std::endl;
}

