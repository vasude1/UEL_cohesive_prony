/* HEADER FOR USERSUBROUTINE UTILITY CLASS
Created by T.Tiirats
------------------------------------------------------------------------------*/


#ifndef  USUB_UTILITYROUTINES_H
#define  USUB_UTILITYROUTINES_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class
// #include <aba_for_c.h>

using namespace Eigen;

class Utility
{
  // PARAMETERS ////////////////////////////////////////////////////////////////
  int nnode;       // Number of nodes
  int ndofel;      // Number of dofs
  int ndim;        // Problem dimension
  int p;           // Order
  int integ_ord;   // Integration order in tangential direction
  double e;        // Thickness // Change to be general !!
  int NSVAR;

 public:
  MatrixXd Coords; // Nodal Coordinates
  MatrixXd u;      // Displacement
  MatrixXd v;      // Displacement
  MatrixXd a;      // Displacement
  MatrixXd du;     // Displacement increment
  double *pMTRX;   // Pointer to output MATRIX array
  double *pVEC;    // Pointer to output VECTOR array
  float time;
  double* SVAR;
  float dt;
  double* parameters;

  // double time;


  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////

  // Constructors
  Utility(int ,int ,int ,int ); // Called from Abaqus

  // SET functions
  void Read_Input(double* COORDS, double* U,double* ,double* ,double* DU, double*,double, double, double*);   // Read Inputs
  void Read_Matrix(double* MATRIX, MatrixXd& Matrix, int m, int n); // Read Matrix from Fort
  void Read_Vector(double* VECTOR, VectorXd& Vector, int n); // Read Vector from Fort
  void Set_Output(double* MAT, double* VEC);    // Output fields
  void Set_Identity();                          // Set matrix to Identity


  // GET functions
  double Get_e(){return e;};
  int Get_p(){return p;};
  int Get_Int_ord(){return integ_ord;};
  int Get_ndofel(){return ndofel;};

  // Write functions
  void Write_Matrix(MatrixXd _K,int m,int n); // Write ouput Matrix
  void Write_Vector(VectorXd _R, int n);      // Write output Vector (Vec*-1 TEMPORARY)
  void set_statevar(double* _SVAR);
  //----------------------------------------------------------------------------


};
#endif
