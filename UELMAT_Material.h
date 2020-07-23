/* HEADER FOR MATERIAL CLASS
Created by T.Tiirats
------------------------------------------------------------------------------*/


// NOTE: IMPLEMENT CONSTRUCTOR FOR ABAQUS MATERIAL ACCESS !!!!!!


#ifndef  UELMAT_MATERIAL_H
#define  UELMAT_MATERIAL_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

#include "UELMAT_ShapeFunctions.h"

using namespace Eigen;

class Material
{
  // INTERNAL PARAMETERS ///////////////////////////////////////////////////////
  MatrixXd Piola_1;
  MatrixXd Piola_2;
  MatrixXd Elasticity;
  MatrixXd DefGrad;
  MatrixXd Green;
  MatrixXd Green_inverse;


  // float c_10 = 0.39E6;
  // float c_01 = 0.97E6;
  float c_10 = 1.0E6;
  float c_01 = 0.0;

  // Variables used outside the class
 public:
  // MatrixXd D;
  //----------------------------------------------------------------------------

  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////
 public:
  // Constructors
  Material(int _dim); // Constructor for Linear Elastic Material

  void Set_Zero();
  //Evaluating functions
  void Set_Material(const ShapeFunc&, const MatrixXd&, const MatrixXd&);  // Build corresponding space
  void Eval_DefGrad(const ShapeFunc&, const MatrixXd&, const MatrixXd&);
  void Eval_Green_and_Inverse();
  void Eval_Elasticity();
  void Eval_Piola_2();
  void Eval_Piola_1();

  //Access functions
  MatrixXd Get_Piola_1() const;
  MatrixXd Get_Piola_2() const;
  MatrixXd Get_Elasticity() const;
  MatrixXd Get_DefGrad() const;


  //----------------------------------------------------------------------------


};
#endif
