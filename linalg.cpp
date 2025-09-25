#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "precice/precice.hpp"
#include <json/json.h>
#include <assert.h>

using namespace precice;

//=========================================================
// Define basic Vector and Matrix objects
//=========================================================
typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Vector&v)
{
   v.resize(dim);
   for(double &d: v)
      d = 0;
}

//=========================================================
// Matrix resize
//=========================================================
void Resize(int dim, Matrix&A)
{
   A.resize(dim);
   for(Vector &vec: A)
      Resize(dim, vec);
}

//=========================================================
// Matrix transpose
//=========================================================
Matrix Transpose(Matrix&A)
{
   int i,j;
   int n = A.size();
   Matrix At;
   Resize(n,At);
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
          At[j][i]= A[i][j];
   }
   return At;
}

//=========================================================
// Matrix-Vector product
//=========================================================
Vector operator*(const Matrix &A, const Vector &x)
{
   int i,j;
   int m = A.size();
   int n = x.size();

   Vector prod(m);

   for(i = 0; i < m; i++)
   {
      prod[i] = 0.;
      for(j = 0; j < n; j++)
          prod[i] += A[i][j]*x[j];
   }
   return prod;
}

//=========================================================
// Matrix-Matrix product
//=========================================================
Matrix operator*(const Matrix &A, const Matrix &B)
{
   assert(A.size() == B.size());  // Assume matrices are square
   int n = A.size();

   int i,j,k;

   Matrix prod;
   Resize(n,prod);

   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         prod[i][j] = 0.;
         for(k = 0; k< n; k++)
            prod[i][j] += A[i][k]*B[k][j];
      }
   }
   return prod;
}

//=========================================================
// Matrix-Matrix-Matrix product for tensor rotation
//=========================================================
void RIRT(const Matrix &R, const Matrix &I0, Matrix &I)
{
   assert(R.size() == I0.size());  // Assume matrices are square
   assert(R.size() == I.size());   // Assume matrices are square
   int n = R.size();
   int i,j,k,l;
   double tmp;
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         I[i][j] = 0.0;
         for(k = 0; k< n; k++)
         {
            tmp = 0.0;
            for(l = 0; l < n; l++)
               tmp += R[j][l]*I0[k][l];
            I[i][j] += R[i][k]*tmp;
         }
      }
   }
}

//=========================================================
// Matrix and Vector print routines
//=========================================================
void print(std::ostream &out, Vector &vec)
{
   out<<"[";
   for(double d: vec)
      out<<d<<" ";
   out<<"]\n";
}

void print(std::ostream &out, Matrix &mat)
{
   out<<"\n";
   for(Vector &vec: mat)
      print(out, vec);
}

//=========================================================
// Matrix and Vector Json conversion routines
//=========================================================
void json2vector(Json::Value& jv, Vector &v)
{
   int dim = jv.size();
   v.resize(dim);
   for (int i = 0; i < dim; ++i)
      v[i] = jv[i].asDouble();
}

void json2matrix(Json::Value& jv, Matrix &m)
{
   int dim = jv.size();
   m.resize(dim);
   for (int i = 0; i < dim; ++i)
   {
      json2vector(jv[i],m[i]);
      assert(m[i].size() == dim);   // Force the matrix to be square
   }
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
         assert(angles_deg.size() == dim);
      }
      else if (ic.isMember("radians"))
      {
         Vector radians;
         json2vector(ic["radians"], radians);
         assert(radians.size() == dim);
      }
      else if (ic.isMember("R"))
      {
         json2matrix(ic["R"], R);
         assert(R.size() == dim);
         // Compute inertia in current configuration
         Matrix id = R*Transpose(R);
         print(std::cout, id);
      }
      else
      {
         Resize(dim, R);
         for (int i = 0;Vector &vec: R)
             vec[i++] = 1.0;
      }

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
      out<<"I0 = ";print(out, I0);

      // State data
      out<<"x = ";print(out, x);
      out<<"v = ";print(out, v);

      out<<"R = ";print(out, R);
      out<<"w = ";print(out, w);

      out<<"I = ";print(out, I);
      Matrix id = R*Transpose(R);
      out<<"RR^t = ";print(std::cout, id);
   }
};

void test(char* rbName)
{
   RigidBody rb(rbName);
   rb.PrintConfig(std::cout);
}

//=========================================================
// The main program
//=========================================================
int main(int argc, char **argv)
{
   // Test -- using default rigid body
   if (argc == 1)
   {
      test(const_cast<char*>("rb.json"));
      return 0;
   }

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

   std::vector<double> boundingBox{
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
       if(participant.requiresWritingCheckpoint())
       {
          rb.saveOldState();
       }
       double preciceDt = participant.getMaxTimeStepSize();
       double solverDt = rb.beginTimeStep();
       double dt = std::min(preciceDt, solverDt);

       participant.readData(meshName, "Displacements", vertexIDs, dt, forces);

       rb.setForces(forces);
       rb.solveTimeStep(dt);
       rb.computeMotion();

       if (displDim != -1)
       {
          participant.writeData(meshName, "Displacements", vertexIDs, rb.computeDisplacements());
       }
       if (velDim != -1)
       {
           participant.writeData(meshName, "Velocities", vertexIDs, rb.computeVelocities());
       }

       participant.advance(dt);

       if (participant.requiresReadingCheckpoint())
       {
          rb.reloadOldState();
       }
       else
       {
          rb.endTimeStep();
       }
   }
   participant.finalize();

   return 0;
}
