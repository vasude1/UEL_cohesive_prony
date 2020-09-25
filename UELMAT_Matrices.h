/* HEADER FOR Matrices CLASS
Created by
------------------------------------------------------------------------------*/


#ifndef  UELMAT_MATRICES_H
#define  UELMAT_MATRICES_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

#include "UELMAT_Material.h"
#include "UELMAT_ShapeFunctions.h"
#include "UELMAT_Integration.h"

using namespace Eigen;
using namespace std;

class Matrices
{
  // INTERNAL PARAMETERS //////////////////////////////////////////////////////
  MatrixXd LHS; // Gausspoint Stiffness contribution
  VectorXd Rhs; // Gausspoint Residual contribution
  MatrixXd Disp_stiff; // Gausspoint Stiffness contribution
  MatrixXd Velo_stiff; // Gausspoint Stiffness contribution
  MatrixXd B_matrix;
  int nnodes;
  double damage;
  double damage_init;
  double damage_fini;
  double separation;
  double cohesive_stiff;
  int number_force;
  int number_elements;
  int number_h;
  int var_per_gp;
  MatrixXd cohesive_prony;
  double velocity_jump;
  Vector2d w;
  double dt_gdt;


//----------------------------------------------------------------------------

  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////
 public:
  // Constructors
  Matrices(int); // Initialize K and R

  void Set_Matrices(const int,const ShapeFunc&,const MatrixXd&,const MatrixXd& ,const MatrixXd&,const MatrixXd&,const MatrixXd&,const MatrixXd&,double* ,float , float , const int* , int , double*);

  void set_Zero();

  void Eval_separation(const VectorXd& ,const MatrixXd& , const int*);

  void Eval_damage(const int,float ,double* , int);

  void Eval_velocity_jump(const MatrixXd&, const MatrixXd&);

  void Eval_stiff(float,const int*,const MatrixXd&,float);

  void Eval_LHS(float, double*);

  void Eval_Rhs(const MatrixXd&,const MatrixXd&,const MatrixXd&,const MatrixXd&,const MatrixXd&,float,float,const VectorXd& ,
              const int* , double*, double*, int, const int);

  void evaluate_BMatrix(MatrixXd& ,const MatrixXd& ,const int* );

  void define_damage_var(float,const MatrixXd& );

  MatrixXd Get_LHS() const;
  MatrixXd Get_Rhs() const;

};
#endif
