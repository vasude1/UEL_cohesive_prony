/* HEADER FOR ASSEMBLY CLASS
Created by T.Tiirats
------------------------------------------------------------------------------*/


#ifndef  UELMAT_ASSEMBLY_H
#define  UELMAT_ASSEMBLY_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

#include "UELMAT_Material.h"
#include "UELMAT_Matrices.h"
#include "UELMAT_ShapeFunctions.h"
#include "UELMAT_Integration.h"

using namespace Eigen;
using namespace std;

class Assembly
{
  // INTERNAL PARAMETERS ///////////////////////////////////////////////////////
  MatrixXd Tangent; // Gausspoint Stiffness contribution
  VectorXd Res; // Gausspoint Residual contribution

//----------------------------------------------------------------------------

  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////
 public:
  // Constructors
  Assembly(int ndofel); // Initialize K and R
  void Eval_Assembly(const MatrixXd&, const VectorXd&, float, float);
  MatrixXd Get_Tangent() const;
  VectorXd Get_Res() const;

};
#endif
