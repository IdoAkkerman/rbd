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
void test(char* fileName)
{
   RigidBody rb(fileName, "rigid_body");

   double dt = rb.GetTimeStep();
   double t = 0;
   Vector forces;
   Resize(2*rb.GetDimension(), forces);
   forces[1] = 1.0;
   forces[2] = 2.0;
   forces[3] = 1.0;
   for (int i = 0; i <100; i++)
   {
      std::cout<<"========================"<<std::endl;
      std::cout<<i<<": t = "<<t<<std::endl;
      std::cout<<"========================"<<std::endl;
      rb.ComputeMotion(dt, forces);
      rb.GetNewState().CheckRotation(1e-8);

      rb.GetNewState().Print(std::cout);
      rb.EndTimeStep();
      t += dt;
   }

   State ref(rb.GetDimension());
   ref.Read(fileName, "ref_condition");
   assert(rb.GetNewState().Equal(ref, 1e-6));
}

Vector computeDisplacements() {Vector disp; return disp;};
Vector computeVelocities() {Vector vel; return vel;};

//=========================================================
// The main program
//=========================================================
int main(int argc, char **argv)
{
   // Test -- using specified rigid body
   if (argc == 2)
   {
      test(argv[1]);
      return 0;
   }

   // Read input
   if (argc != 5)
   {
      std::cout<<"Wrong number of inputs.\n";
      std::cout<<"  rbd rbFile\n";
      std::cout<<"  rbd rbFile configFile solverName meshName\n";
      return -1;
   }

   std::string rbFile(argv[1]);
   std::string configFile(argv[2]);
   std::string solverName(argv[3]);
   std::string meshName(argv[4]);

   // Configure rigid body
   RigidBody rb(rbFile);

   // Configure precice
   Participant participant(solverName, configFile, 0, 1);

   std::vector<double> boundingBox
   {
      -9999999, 9999999, // x-axis min and max
         -9999999, 9999999, // y-axis min and max
         -9999999, 9999999  // z-axis min and max
      };

   // Define region of interest, where we want to obtain the direct access.
   // See also the API documentation of this function for further notes.
   participant.setMeshAccessRegion(meshName, boundingBox);

   // initialize preCICE as usual
   participant.initialize();

   const int dimensions = participant.getMeshDimensions(meshName);
   const int vertexSize = participant.getMeshVertexSize(meshName);

   std::vector<double>  coordinates(vertexSize * dimensions);
   std::vector<int>     vertexIDs(vertexSize);

   participant.getMeshVertexIDsAndCoordinates(meshName,
                                              vertexIDs,
                                              coordinates);

   const int forceDim = participant.getDataDimensions(meshName, "Forces");
   const int displDim = participant.getDataDimensions(meshName, "Displacements");
   const int velDim = participant.getDataDimensions(meshName, "Velocities");

   std::vector<double> forces(vertexSize*dimensions);

   // Main timestepping/iteration loop
   while (participant.isCouplingOngoing())
   {
      if (participant.requiresWritingCheckpoint())
      {
         rb.SaveOldState();
      }
      double preciceDt = participant.getMaxTimeStepSize();
      double solverDt = rb.GetTimeStep();
      double dt = std::min(preciceDt, solverDt);

      participant.readData(meshName, "Displacements", vertexIDs, dt, forces);

      rb.ComputeMotion(dt, forces);

      if (displDim != -1)
      {
         participant.writeData(meshName, "Displacements", vertexIDs,
                               computeDisplacements());
      }
      if (velDim != -1)
      {
         participant.writeData(meshName, "Velocities", vertexIDs, computeVelocities());
      }

      participant.advance(dt);

      if (participant.requiresReadingCheckpoint())
      {
         rb.ReloadOldState();
      }
      else
      {
         rb.EndTimeStep();
      }
   }
   participant.finalize();

   return 0;
}
