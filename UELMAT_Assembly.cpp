/* ASSEMBLY CLASS TO DEFINE UELMAT ELMENT LHS AND RHS
Created by T.Tiirats
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Assembly.h"
#include "UELMAT_Material.h"
#include "UELMAT_Matrices.h"
#include "UELMAT_ShapeFunctions.h"
#include "UELMAT_Integration.h"

// Namespaces
using namespace Eigen;
using namespace std;



// CLASS MATERIAL FUNCTIONS ///////////////////////////////////////////////////


// Constructors ****************************************************************
Assembly::Assembly(int ndofel){

  Tangent.resize(ndofel,ndofel);     // Stiffness matrix
  Tangent.setZero();
  Res.resize(ndofel);            // RHS
  Res.setZero();
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



// Eval functions **************************************************************
void Assembly::Eval_Assembly(const MatrixXd& LHS, const VectorXd& Res_, float weight, float jac_det){

  Tangent += weight*jac_det*LHS;
  // cout<<Tangent<<endl;

  Res += weight*jac_det*Res_;
  // cout<<Res<<endl;
};
//------------------------------------------------------------------------------

MatrixXd Assembly::Get_Tangent() const{
  return Tangent;

};

VectorXd Assembly::Get_Res() const{
  return Res;

};
