/* HEADER FOR SHAPEFUNCTION CLASS
Created by T.Tiirats
------------------------------------------------------------------------------*/

#ifndef  UELMAT_ShapeFunctions_H
#define  UELMAT_ShapeFunctions_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

#include "UELMAT_Integration.h"// Eigen class

using namespace Eigen;

class ShapeFunc
{
  // INTERNAL PARAMETERS ///////////////////////////////////////////////////////
  int n_shape=2; // Number of shape fucntions, 3 for tri, 4 for quad in 2d
  int dim=1; // Dimension of the problem, 2d for now
  Vector2d Shape; //Shape functions
  Vector2d DShape_parent; // Derivatives in parent space
  // MatrixXd DShape_physical; // Derivatives in physical space
  double Jacobian; //Jacbian of transformation between parent andphysical spaces
  Vector2d tangent_deformed;  //Tangent vector in the deformed configuration
  Vector2d normal_deformed;   //Normal vector in the deformed configuration


  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////
 public:
  // Constructors
  ShapeFunc(int _n_shape, int _dim); //Change to appropriate shape

  void Set_Zero();
  // Set functions
  void Set_Approx(const float, const MatrixXd&, const MatrixXd&, const int*);  // Build corresponding space

  // Eval functions
  void Eval_Shape(const float);
  //-----------------------------------------------------------------------------
  void Eval_DShape_parent();
  void Eval_Jacobian(const MatrixXd&);
  void Eval_tangent(const MatrixXd&, const MatrixXd&, const int*);
  void Eval_normal(const Vector2d&);


  //Access functions
  Vector2d Get_Shape() const;
  float Get_Jacobian_det() const;
  Vector2d Get_tangent() const;
  Vector2d Get_normal() const;

};

#endif
