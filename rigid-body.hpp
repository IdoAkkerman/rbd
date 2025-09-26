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
class State
{
public:
   int dim;
   Vector x,v;  // Displacement & Velocity;
   Matrix R;    // Rotation matrix
   Vector w;    // Angular velocity

   Matrix I;    // Inertia tensor --> in current configuration

   State(){};

   State(State &org)
   {
      dim = org.dim;
      x = org.x;
      v = org.v;
      R = org.R;
      w = org.w;
      I = org.I;
   }

   void Initialize(int dim_)
   {
       dim = dim_;
       Resize(dim, x);
       Resize(dim, v);
       Resize(dim, w);
       Resize(dim, R);
       for (int i = 0;Vector &vec: R)
          vec[i++] = 1.0;
       Resize(dim, I);
   }

   void SetInertiaTensor(Matrix &I0)
   {
      assert(I0.size() == dim);
      I = R*I0*Transpose(R);
   }

   void CheckRotation()
   {
      // Check rotation matrix
      Matrix id = R*Transpose(R);

      for (int i = 0; Vector vec: id)
         vec[i++] -= 1.0;
      double norm = Norm(id);
      if (norm > 1e-10)
      {
         cout<<" Norm = " << norm<<std::endl
         Matrix id = R*Transpose(R);
         ::Print(std::cout, id);
      }
      //assert(norm < 1e-12);
   }

   void Read(Json::Value &data)
   {
      // Read state data -- zero if not defined
      // Displacement
      if (data.isMember("x"))
      {
         json2vector(data["x"], x);
         assert(x.size() == dim);
      }
      else
      {
         Resize(dim, x);
      }

      // Velocity
      if (data.isMember("v"))
      {
         json2vector(data["v"], v);
         assert(v.size() == dim);
      }
      else
      {
         Resize(dim, v);
      }

      // Angular Velocity
      if (data.isMember("w"))
      {
         json2vector(data["w"], w);
         assert(w.size() == (dim*(dim-1))/2);
      }
      else
      {
         Resize((dim*(dim-1))/2, w);
      }

      // Rotation
      if (data.isMember("angles_deg"))
      {
         Vector angles_deg;
         json2vector(data["angles_deg"], angles_deg);
         assert(angles_deg.size() == (dim*(dim-1))/2);
         for (double &d: angles_deg)
            d *= M_PI/180;

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
         for (int i = 0;Vector &vec: R)
             vec[i++] = 1.0;
      }
   };

   void Print(std::ostream &out)
   {
      out<<"x = ";::Print(out, x);
      out<<"v = ";::Print(out, v);

      out<<"R = ";::Print(out, R);
      out<<"w = ";::Print(out, w);
   }
};

//=========================================================
// Define the rigid body class
//=========================================================
class RigidBody
{
private:
   int dim;     // Dimension of the problem (2D/3D)

   double m;    // Mass
   Matrix I0;    // Mass moment of inertia tensor

   double dt_max;

   State state_old, state_new;

   void Construct(Json::Value &data)
   {
      // Read object data
      dim = data["dimension"].asInt();
      m = data["mass"].asDouble();
      json2matrix(data["I0"],I0);
      assert(I0.size() == dim);

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

   int GetDim(){return dim;};

   double beginTimeStep(){return dt_max;};

   void computeMotion(double dt, Vector &forces)
   {
      state_old.Print(std::cout);
      state_new.Print(std::cout);
       Vector F(dim);
       Vector M(dim);
       int i = 0;
       for(double &d: F)
          d = forces[i++];
       for(double &d: M)
          d = forces[i++];

       state_new.v =  state_old.v + (dt/m)*F;
       state_new.x =  state_old.x + (dt/2)*(state_new.v+state_old.v);

       Matrix eye;
       Resize(dim, eye);
       for(int i = 0; Vector & vec: eye)
          vec[i++] = 1.0;

       for(int i = 0; i < 112; i++)
       {
          std::cout<<"i = "<<i<<std::endl;
          state_new.SetInertiaTensor(I0);
          Matrix I_inv = Inverse(state_new.I);
          state_new.w =  I_inv*(state_old.I*state_old.w + dt*I_inv*M);

          Matrix Q = Skew(0.5*(state_new.w+state_old.w));
          Matrix Qp = eye + (dt/2)*Q;
          Matrix Qm_inv = Inverse(eye + (-dt/2)*Q);

          Matrix Rn= Qm_inv*Qp*state_old.R;
          double norm = Norm(state_new.R +(-1.0)*Rn);
          state_new.R =  Rn;
          ::Print(std::cout, state_new.R );

          if(norm < 1e-8) break;
       }
       state_old.Print(std::cout);
       state_new.Print(std::cout);
   };

   void reloadOldState(){}; // Not needed
   void saveOldState(){};  // Not needed
   void endTimeStep(){ state_old = state_new; };

   State& GetNewState(){ return state_new; };
   State& GetOldState(){ return state_old; };
   
   void Print(std::ostream &out)
   {
      // Object data
      out<<"Dimension = "<<dim<<std::endl;
      out<<"Mass = "<<m<<std::endl;
      out<<"I0 = ";::Print(out, I0);

      // State data
      state_old.Print(out);
      state_new.Print(out);
   }
};

#endif
