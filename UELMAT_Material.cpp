/* MATERIAL CLASS TO DEFINE UELMAT ELMENT INTEGRATION
Created by T.Tiirats
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Material.h"

// Namespaces
using namespace Eigen;
using namespace std;


// Constructors ****************************************************************
Material::Material(int _dim){

  Piola_1.resize(_dim,_dim);
  Piola_2.resize(_dim,_dim);

  DefGrad.resize(_dim,_dim);
  Green.resize(_dim,_dim);
  Green_inverse.resize(_dim,_dim);

  Elasticity.resize(3,3);
};

void Material::Set_Material(const ShapeFunc& Shape, const MatrixXd& _U, const MatrixXd& _Coords){
  //Do not change the order
  Set_Zero();
  //Evaluate deformation gradient
  Eval_DefGrad(Shape, _U, _Coords);

  //Evaluate Green strain using DefGrad above
  Eval_Green_and_Inverse();

  //Computes Elasticity tensor using Green strain and inverse above
  Eval_Elasticity();

  //Compute PK2 - uses Green strain and inverse
  Eval_Piola_2();

  //Compute PK1 - uses PK2
  Eval_Piola_1();
};

void Material::Set_Zero(){
  Piola_1.setZero();
  Piola_2.setZero();

  DefGrad.setZero();
  Green.setZero();
  Green_inverse.setZero();

  Elasticity.setZero();

};
//------------------------------------------------------------------------------
// Evaluate Deformation gradient
void Material::Eval_DefGrad(const ShapeFunc& Shape, const MatrixXd& _U, const MatrixXd& _Coords){
  MatrixXd new_position = _U + _Coords;
  DefGrad = new_position.transpose() * Shape.Get_DShape_physical();
};

// Evaluate Deformation gradient
void Material::Eval_Green_and_Inverse(){
  Green = DefGrad.transpose()*DefGrad;
  Green_inverse = Green.inverse();
};



// Evaluate Elasticity tensor 2d
void Material::Eval_Elasticity(){

  float I3 = Green.determinant();
  float I1 = Green.trace();
  float scalar = 4.0*(c_10/I3 + c_01/I3*I1);

  Elasticity(0,0) = 1.0-1.0;
  Elasticity(0,1) = 1.0;
  Elasticity(0,2) = 0;
  Elasticity(1,0) = 1.0;
  Elasticity(1,1) = 1.0-1.0;
  Elasticity(1,2) = 0;
  Elasticity(2,0) = 0;
  Elasticity(2,1) = 0;
  Elasticity(2,2) = 0.0-1.0/2.0;
  Elasticity = 4.0*c_01*Elasticity;

  Elasticity(0,0) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(0,0);
  Elasticity(0,1) += 0.5*scalar*Green_inverse(0,1)*Green_inverse(0,1);
  Elasticity(0,2) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(0,1);
  Elasticity(1,0) += 0.5*scalar*Green_inverse(1,0)*Green_inverse(1,0);
  Elasticity(1,1) += 0.5*scalar*Green_inverse(1,1)*Green_inverse(1,1);
  Elasticity(1,2) += 0.5*scalar*Green_inverse(1,0)*Green_inverse(1,1);
  Elasticity(2,0) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(1,0);
  Elasticity(2,1) += 0.5*scalar*Green_inverse(0,1)*Green_inverse(1,1);
  Elasticity(2,2) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(1,1);

  Elasticity(0,0) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(0,0);
  Elasticity(0,1) += 0.5*scalar*Green_inverse(0,1)*Green_inverse(0,1);
  Elasticity(0,2) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(0,1);
  Elasticity(1,0) += 0.5*scalar*Green_inverse(1,0)*Green_inverse(1,0);
  Elasticity(1,1) += 0.5*scalar*Green_inverse(1,1)*Green_inverse(1,1);
  Elasticity(1,2) += 0.5*scalar*Green_inverse(1,0)*Green_inverse(1,1);
  Elasticity(2,0) += 0.5*scalar*Green_inverse(0,0)*Green_inverse(1,0);
  Elasticity(2,1) += 0.5*scalar*Green_inverse(0,1)*Green_inverse(1,1);
  Elasticity(2,2) += 0.5*scalar*Green_inverse(0,1)*Green_inverse(1,0);


  Elasticity(0,0) += scalar*Green_inverse(0,0)*Green_inverse(0,0);
  Elasticity(0,1) += scalar*Green_inverse(0,0)*Green_inverse(1,1);
  Elasticity(0,2) += scalar*Green_inverse(0,0)*Green_inverse(0,1);
  Elasticity(1,0) += scalar*Green_inverse(1,1)*Green_inverse(0,0);
  Elasticity(1,1) += scalar*Green_inverse(1,1)*Green_inverse(1,1);
  Elasticity(1,2) += scalar*Green_inverse(1,1)*Green_inverse(0,1);
  Elasticity(2,0) += scalar*Green_inverse(0,1)*Green_inverse(0,0);
  Elasticity(2,1) += scalar*Green_inverse(0,1)*Green_inverse(1,1);
  Elasticity(2,2) += scalar*Green_inverse(0,1)*Green_inverse(0,1);


  Elasticity(0,0) -= 4.0*c_01/I3*Green_inverse(0,0);
  Elasticity(0,1) -= 4.0*c_01/I3*Green_inverse(0,0);

  Elasticity(1,0) -= 4.0*c_01/I3*Green_inverse(1,1);
  Elasticity(1,1) -= 4.0*c_01/I3*Green_inverse(1,1);

  Elasticity(2,0) -= 4.0*c_01/I3*Green_inverse(0,1);
  Elasticity(2,1) -= 4.0*c_01/I3*Green_inverse(0,1);

  Elasticity(0,0) -= 4.0*c_01/I3*Green_inverse(0,0);
  Elasticity(0,1) -= 4.0*c_01/I3*Green_inverse(1,1);
  Elasticity(0,2) -= 4.0*c_01/I3*Green_inverse(0,1);

  Elasticity(1,0) -= 4.0*c_01/I3*Green_inverse(0,0);
  Elasticity(1,1) -= 4.0*c_01/I3*Green_inverse(1,1);
  Elasticity(1,2) -= 4.0*c_01/I3*Green_inverse(0,1);

};

//Evaluate the PK2 stress
void Material::Eval_Piola_2(){

  float I3 = Green.determinant();
  float I1 = Green.trace();
  Matrix2d identity = Matrix2d::Identity();
  float scalar = 2*(c_10/I3+c_01*I1/I3);
  Piola_2 = 2*c_10*identity;
  Piola_2 += 2*c_01*((I1+1.0/I3)*identity-Green);
  Piola_2 -= scalar*Green_inverse;

};

//Evaluate the PK1 stress
void Material::Eval_Piola_1(){

  Piola_1 = DefGrad*Piola_2;
};


//******************************************************************************
//Access functions
//Returns PK1
MatrixXd Material::Get_Piola_1() const{

  return Piola_1;
};

//Returns PK2
MatrixXd Material::Get_Piola_2() const{

  return Piola_2;
};


//Returns Elasticity tensor
MatrixXd Material::Get_Elasticity() const{

  return Elasticity;
};

MatrixXd Material::Get_DefGrad() const{

  return DefGrad;
};
