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
#include "linalg.hpp"

//=========================================================
// Convert thte angles to a rotation matrix
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
class RigidBody
{
private:
   int dim;     // Dimension of the problem (2D/3D)

   double m;    // Mass
   Matrix I0;   // Mass moment of inertia tensor
   Matrix I;    // Mass moment of inertia tensor

   Vector x,v;  // Displacement & Velocity;
   Matrix R;    // Rotation matrix
   Vector w;    // Angular velocity

   Vector disp, vel;

   void Construct(Json::Value &data)
   {
      // Read object data
      dim = data["dimension"].asInt();
      m = data["mass"].asDouble();
      json2matrix(data["I0"],I0);
      assert(I0.size() == dim);

      Json::Value ic = data["initial_condition"];
 
      // Read state data -- zero if not defined
      // Displacement
      if (ic.isMember("x"))
      {
         json2vector(ic["x"], x);
         assert(x.size() == dim);
      }
      else
      {
         Resize(dim, x);
      }

      // Velocity
      if (ic.isMember("v"))
      {
         json2vector(ic["v"], v);
         assert(v.size() == dim);
      }
      else
      {
         Resize(dim, v);
      }

      // Angular Velocity
      if (ic.isMember("w"))
      {
         json2vector(ic["w"], w);
         assert(w.size() == dim);
      }
      else
      {
         Resize(dim, w);
      }

      // Rotation
      if (ic.isMember("angles_deg"))
      {
         Vector angles_deg;
         json2vector(ic["angles_deg"], angles_deg);
         assert(angles_deg.size() == (dim*(dim-1))/2);
         for (double &d: angles_deg)
            d *= M_PI/180;

         R = GetRotation(angles_deg);
      }
      else if (ic.isMember("radians"))
      {
         Vector radians;
         json2vector(ic["radians"], radians);
         assert(radians.size()  == (dim*(dim-1))/2);
         R = GetRotation(radians);
      }
      else if (ic.isMember("R"))
      {
         json2matrix(ic["R"], R);
         assert(R.size() == dim);
      }
      else
      {
         Resize(dim, R);
         for (int i = 0;Vector &vec: R)
             vec[i++] = 1.0;
      }

      // Check rotation matrix
      Matrix id = R*Transpose(R);

      double frobenius = 0.0;
      for (int i = 0; Vector vec: id)
      {
         vec[i++] -= 1.0;
         for (double d: vec)
            frobenius += d*d;
      }
      if (frobenius  > 1e-12) Print(std::cout, id);
      assert(frobenius  < 1e-12);
      
      // Compute inertia in current configuration
      Resize(dim, I);
      RIRT(R,I0,I);
   }

public:

   RigidBody(std::string fileName, std::string rbName)
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

   RigidBody(std::string fileName)
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

   RigidBody(Json::Value &data)
   {
      Construct(data);
   }

   void setForces(Vector forces){};
   double beginTimeStep(){return 0.0;};
   void solveTimeStep(double dt){};

   void computeMotion(){};
   Vector computeDisplacements(){return disp;};
   Vector computeVelocities(){return vel;};

   void reloadOldState(){};
   void endTimeStep(){};
   void saveOldState(){};

   void PrintConfig(std::ostream &out)
   {
      // Object data
      out<<"Dimension = "<<dim<<std::endl;
      out<<"Mass = "<<m<<std::endl;
      out<<"I0 = ";Print(out, I0);

      // State data
      out<<"x = ";Print(out, x);
      out<<"v = ";Print(out, v);

      out<<"R = ";Print(out, R);
      out<<"w = ";Print(out, w);

      out<<"I = ";Print(out, I);
      Matrix id = R*Transpose(R);
      out<<"RR^t = ";Print(std::cout, id);
   }
};

#endif
