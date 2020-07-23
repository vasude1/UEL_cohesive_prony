/* SHAPEFUNCTION CLASS TO DEFINE UELMAT ELMENT INTEGRATION
Created by VaKa
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
//#include "UELMAT_Integration.h"
#include "UELMAT_ShapeFunctions.h"

// Namespaces
using namespace Eigen;
using namespace std;



// CLASS SHAPEFUNC FUNCTIONS ///////////////////////////////////////////////////


// Constructors ****************************************************************
ShapeFunc::ShapeFunc(int _n_shape = 2,int _dim = 1){
  n_shape = _n_shape;
  dim = _dim;

};
//------------------------------------------------------------------------------

// Set functions ***************************************************************

// Build correct space
void ShapeFunc::Set_Approx(const float _Gp, const MatrixXd& _Coords, const MatrixXd& _U, const int* sequence){
  // Build lagrangian space:
  //cout << " Using Linear Lagrangian Approximation Space. " << endl;
  Set_Zero();
  // std::cout << _Coords<< '\n';
  // Construct Shapefunctions vector
  Eval_Shape(_Gp);

  //Evaluate the Jacobian
  Eval_Jacobian(_Coords);

  Eval_tangent(_Coords, _U, sequence);

  Eval_normal(tangent_deformed);

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void ShapeFunc::Set_Zero(){
  Shape.setZero(); //Shape functions
  DShape_parent.setZero(); // Derivatives in parent space
  tangent_deformed.setZero();
  normal_deformed.setZero();
  Jacobian = 0.0;

};

// Eval functions **************************************************************

//Evaluate shape functions
void ShapeFunc::Eval_Shape(const float _Gp){
  float Gp_e = _Gp;
  Shape(0) = (1-Gp_e)/2.0;
  Shape(1) = (1+Gp_e)/2.0;
};


//Evaluate derivatives in parent element
void ShapeFunc::Eval_DShape_parent(){
  DShape_parent(0) = -1.0/2.0;
  DShape_parent(1) = 1.0/2.0;
 };


// Eval Jacobian
void ShapeFunc::Eval_Jacobian(const MatrixXd& _Coords){

  Jacobian = 0.5*abs(_Coords(1,0)-_Coords(0,0));

};

void ShapeFunc::Eval_tangent(const MatrixXd& _Coords, const MatrixXd& U, const int* sequence)
{
  MatrixXd deformed_coords = _Coords+U;
  MatrixXd left_end = 0.5*(deformed_coords(*(sequence+3),seq(0,1))+deformed_coords(*(sequence+0),seq(0,1)));
  MatrixXd right_end = 0.5*(deformed_coords(*(sequence+2),seq(0,1))+deformed_coords(*(sequence+1),seq(0,1)));

  tangent_deformed = (right_end-left_end).transpose();
  tangent_deformed /= tangent_deformed.norm();

};

void ShapeFunc::Eval_normal(const Vector2d& tangent_deformed)
{
  normal_deformed(0) = -tangent_deformed(1);
  normal_deformed(1) = tangent_deformed(0);

};


//Access functions************************************************************
//Gives shape functions
Vector2d ShapeFunc::Get_Shape() const{
  return Shape;
};

//Gives the determinant of Jacobian
float ShapeFunc::Get_Jacobian_det() const{
  return Jacobian;
};

Vector2d ShapeFunc::Get_tangent() const{
  return tangent_deformed;
};

Vector2d ShapeFunc::Get_normal() const{
  return normal_deformed;
};
