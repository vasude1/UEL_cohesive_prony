/* USERSUBROUTINES FOR RTLM
Created by T.Tiirats
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <Eigen/Dense> // Eigen class
// #include <aba_for_c.h>

// OWN INCLUDES
#include "USUB_UtilityRoutines.h"


// Namespaces
using namespace Eigen;
using namespace std;


// CLASS UTILITY FUNCTIONS /////////////////////////////////////////////////////

// Constructor
Utility::Utility(int NNODE,int NDOFEL,int integ_ord,int NSVARS){
  nnode = NNODE;
  ndofel = NDOFEL;
  ndim = ndofel/nnode;
  p = nnode/2-1;
  // integ_ord = _integ_ord;
  NSVAR = NSVARS;

  Coords.resize(nnode,ndim);
  u.resize(nnode,ndim);
  du.resize(nnode,ndim);
  v.resize(nnode,ndim);
  a.resize(nnode,ndim);

  // SVAR.resize(NSVAR,1);

};
//------------------------------------------------------------------------------


// Read fortran arrays (from column-major to row-major)
void Utility::Read_Input(double* COORDS,double* U,double* V, double* A ,double* DU, double* _SVAR, double DTIME, double TIME,
                              double* PARAMS){
  Read_Matrix(COORDS,Coords,nnode,ndim); // Read Coordinates Matrix
  Read_Matrix(U,u,nnode,ndim);               // Read U_{i-1}
  Read_Matrix(V,v,nnode,ndim);               // Read U_{i-1}
  Read_Matrix(A,a,nnode,ndim);               // Read U_{i-1)
  Read_Matrix(DU,du,nnode,ndim);               // Read U_{i-1}
  // Read_Vector(_SVAR,SVAR,NSVAR);
  SVAR = _SVAR;
  dt = DTIME;
  time = TIME;
  parameters = PARAMS;
};
//------------------------------------------------------------------------------

void Utility::Read_Matrix(double* MATRIX, MatrixXd& Matrix, int m, int n){
  int k=0;
  for (int i=0;i<m;++i){
    for (int j=0;j<n;++j){
      Matrix(i,j) = MATRIX[k];
      k += 1;
    };
  };

};
//------------------------------------------------------------------------------

void Utility::Read_Vector(double* VECTOR, VectorXd& Vector, int n){
  for (int i=0;i<n;++i){
    Vector(i) = VECTOR[i];
    // cout<<Vector(i)<<endl;
  };

};
//------------------------------------------------------------------------------

// Set output fortran array pointer
void Utility::Set_Output(double* MAT,double* VEC){
  pVEC = VEC;
  pMTRX = MAT;
  // std::cout << "SVAR" << '\n';
  // std::cout << *SVARS << '\n';
};
//------------------------------------------------------------------------------

// Set identity
void Utility::Set_Identity(){
  MatrixXd I = MatrixXd::Identity(ndofel, ndofel);

};
//------------------------------------------------------------------------------

// Write fortran array (from row-major to column-major)
void Utility::Write_Matrix(MatrixXd _K,int m,int n){
  int k=0;
  for (int j=0;j<n;++j){
    for (int i=0;i<m;++i){
      pMTRX[k] = _K(i,j);
      k += 1;
    };
  };

};
//------------------------------------------------------------------------------

// Write fortran list (from vector to vector)
void Utility::Write_Vector(VectorXd _R,int n){
  for (int i=0;i<n;++i){
    pVEC[i] =_R(i);
  };

};
//------------------------------------------------------------------------------

void Utility::set_statevar(double* _SVAR){
  *SVAR = *_SVAR;
};
