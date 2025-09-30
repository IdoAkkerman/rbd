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
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "precice/precice.hpp"
#include <json/json.h>
#include <assert.h>

#include "vector.hpp"
#include "matrix.hpp"
#include "state.hpp"
#include "rigid-body.hpp"

using namespace precice;

//=========================================================
// The test program
//=========================================================
int main(int argc, char **argv)
{
   // Test -- using specified rigid body
   if (argc != 2)
   {

      return -3;
   }

   std::string fileName(argv[1]);
   RigidBody rb(fileName, "rigid_body");
   rb.SetVerbose();

   Vector forces;
   Resize(rb.GetDimension(), forces);
   for (int i = 0; i <forces.size(); i++) { forces[i] = i+1; }

   Vector moments;
   Resize(rb.GetRDimension(), moments);
   for (int i = 0; i <moments.size(); i++) { moments[i] = i+1; }

   std::fstream fs;
   fs.open ("output.txt", std::fstream::out);
   rb.GetOldState().PrintHeader(fs);
   rb.GetOldState().Print(0.0,fs);
   rb.GetOldState().PrettyPrint(std::cout);
   for (int i = 0; i <1000; i++)
   {
      double dt = rb.GetTimeStep();
      std::cout<<"==========================================="<<std::endl;
      std::cout<<i+1<<" "<<dt<<": "<<rb.GetOldState().t
               <<" --> "<<rb.GetOldState().t+dt<<std::endl;
      std::cout<<"==========================================="<<std::endl;
      rb.ComputeMotion(dt, forces, moments);
      rb.GetNewState().CheckRotation(1e-12);

      rb.EndTimeStep();

      rb.GetNewState().PrettyPrint(std::cout);
      rb.GetOldState().Print(dt,fs);
   }

   State ref(rb.GetDimension());
   bool hasRef = ref.Read(fileName, "ref_condition");
   if (hasRef)
   {
      // Compare
      bool isEqual = rb.GetNewState().Equal(ref, 1e-12);
      if (isEqual)
      {
         std::cout<<"Solution matches reference!! \n";
      }
      else
      {
         return -1;
      }
   }
   else
   {
      std::cout<<"\"ref_condition\" : \n";
      rb.GetNewState().PrintJson(std::cout);
      std::cout<<std::flush;
      return -2;
   }
   return 0;
}

