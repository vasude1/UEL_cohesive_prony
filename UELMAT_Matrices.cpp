/* Matrices CLASS TO DEFINE UELMAT ELMENT LHS AND RHS
Created by Va Ka
--------------------------------------------------------------------------------
Last edit:
Changed RHS - removed if conditionsl
changed LHS - removed if conditional
changed the stiffness in LHS - added conditions for damage
Using inconsistent stiffness matrix for the velocity dependent terms
Using delta_0= 1E-6, delta_f = 3E-3, eta = 7E6
------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Matrices.h"
#include "UELMAT_Material.h"
#include "UELMAT_ShapeFunctions.h"
#include "UELMAT_Integration.h"

// Namespaces
using namespace Eigen;
using namespace std;

// Signature


// CLASS MATERIAL FUNCTIONS ///////////////////////////////////////////////////


// Constructors ****************************************************************
Matrices::Matrices(int ndofel){

  LHS.resize(ndofel,ndofel); // Gausspoint Stiffness contribution

  Disp_stiff.resize(ndofel,ndofel);

  Velo_stiff.resize(ndofel,ndofel);

  Rhs.resize(ndofel);

  nnodes = ndofel/2;

  cohesive_prony << 0.096092126464814,	1E-14,
  0.131669775836881,	1E-13,
  0.229532673957081,	1E-12,
  0.231886167949743,	1E-11,
  0.183771286292404,	1E-10,
  0.072560462039951,	1E-09,
  0.025498195885329,	1E-08,
  0.009082879762867,	1E-07,
  0.005070356052601,	1E-06,
  0.003003882953384,	1E-05,
  0.002091305180757,	0.0001,
  0.00355114288563,0.001;

  // cohesive_prony << 0.0,	1E-14,
  // 0.0,	1E-13,
  // 0.0,	1E-12,
  // 0.0,	1E-11,
  // 0.0,	1E-10,
  // 0.0,	1E-09,
  // 0.0,	1E-08,
  // 0.0,	1E-07,
  // 0.0,	1E-06,
  // 0.0,	1E-05,
  // 0.0,	0.0001,
  // 0.0,0.001;

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void Matrices::Set_Matrices(const int lflag,const ShapeFunc& Shape,const MatrixXd& Coords,const MatrixXd& U,const MatrixXd& DU,const MatrixXd& V,const MatrixXd& A,const MatrixXd& Mass ,double* SVARS,float delta_t,
                          float time, const int* sequence, int gp, double* hht)
{
  MatrixXd shape_funcs = Shape.Get_Shape();
  MatrixXd tangent_deformed = Shape.Get_tangent();
  MatrixXd normal_deformed = Shape.Get_normal();

  cohesive_stiff = 1E13;

  evaluate_BMatrix(B_matrix, shape_funcs, sequence);

  define_damage_var(time,Coords);

  set_Zero();

  Eval_separation(shape_funcs,U,normal_deformed,sequence);

  Eval_damage(lflag,separation,SVARS,gp);

  Eval_stiff(shape_funcs,tangent_deformed,normal_deformed, separation, sequence, U,V, delta_t);

  Eval_LHS(time,hht);

  Eval_Rhs(Mass,U,DU,V,A,time,delta_t, shape_funcs, sequence, hht,SVARS,gp,lflag);

  // throw;
};

void Matrices::set_Zero(){
  LHS.setZero();
  Rhs.setZero();
  Disp_stiff.setZero();
  Velo_stiff.setZero();
};

void Matrices::define_damage_var(float time,const MatrixXd& Coords ){

  // Assign the actual damage parameters and initialize the seed crack
  float x_centroid = real(Coords(0,0)+Coords(1,0)+Coords(2,0)+Coords(3,0));
  x_centroid /= 4.0;
  damage_init = 1E-6; // 1E-6
  damage_fini = 2E-3; //2E-3
  if (time>1.0)   // time>1.0
  {

    if(x_centroid < 1.5E-3)  // Seed crack
    {
      damage_init = 1E-7;
      damage_fini = 3E-3;
    }
  }

};

void Matrices::Eval_velocity_jump(const MatrixXd& B_matrix, const MatrixXd& V){
  VectorXd V_col= V.reshaped<RowMajor>();
  VectorXd jump = B_matrix.transpose()*V_col;
  double jump_x = real(jump(0));
  double jump_y = real(jump(1));
  velocity_jump = pow(pow(jump_x,2)+pow(jump_y,2),0.5);

};

void Matrices::Eval_separation(const VectorXd& shape_funcs,const MatrixXd& U,const Vector2d& normal_deformed, const int* sequence){

   MatrixXd left_separation;  // Separation on the Left
  left_separation = U(*(sequence+3),seq(0,1)) - U(*(sequence+0),seq(0,1));

   MatrixXd right_separation; //Separation on the right end
  right_separation = U(*(sequence+2),seq(0,1)) - U(*(sequence+1),seq(0,1));

   MatrixXd gp_separation;  // Separation at Gauss point = N1*u1+N2*u2
  gp_separation = shape_funcs(0)*left_separation+shape_funcs(1)*right_separation;

  float separation_x = real(gp_separation(0,0));
  float separation_y = real(gp_separation(0,1));
  separation = pow(pow(separation_x,2)+pow(separation_y,2),0.5);

};


void Matrices::Eval_damage(const int lflag, float separation,double* SVARS, int gp){

  // If fully damaged, return
  damage=0.0;
  if(separation > damage_init)
  {
    damage = damage_fini*(separation-damage_init)/(separation)/(damage_fini-damage_init);
  }

  if(damage > *(SVARS+gp))
  {
    if(damage>=1)
    {
      damage=1.0;
    }
    if(lflag==1)
    {
        *(SVARS+gp)=damage;
    }
    return;
  }

  if(damage < *(SVARS+gp))
  {
    if(lflag==1)
    {
        damage=*(SVARS+gp);
    }
    return;
  }
  return;

};


void Matrices::Eval_stiff(const VectorXd& shape_funcs,const Vector2d& tangent_deformed,const Vector2d& normal_deformed,
                  float separation,const int* sequence, const MatrixXd& U, const MatrixXd& V, float delta_t){

    // Elemental matrices
    MatrixXd stiff(2,2);
    stiff.setZero();

    // To determine if the faces are opening or closing
    VectorXd V_col= V.reshaped<RowMajor>();
    double multiplication_factor=1;

    for(int i=0;i<number_elements;++i)
    {
      multiplication_factor += cohesive_prony(i,0)*exp(-delta_t/2.0/cohesive_prony(i,1))/(1.0-(cohesive_prony.colwise().sum())(0));
      // std::cout << "Begin Print" << '\n';
      // std::cout << delta_t << '\n';
      // std::cout << cohesive_prony(i,1) << '\n';
      // std::cout << cohesive_prony(i,0)*exp(-delta_t/2.0/cohesive_prony(i,1)) << '\n';
    }
    // std::cout << multiplication_factor << '\n';
    // throw;
    // std::cout << exp(-delta_t/2.0*1000.0) << '\n';
    stiff(0,0) = cohesive_stiff*multiplication_factor;
    stiff(1,1) = cohesive_stiff*multiplication_factor;

    stiff(0,0) *= (1-damage);
    stiff(1,1) *= (1-damage);

    // if(damage>0)
    // {
    //   if((B_matrix.transpose()*V_col)[1]<0.0) // Closing faces
    //   {
    //     stiff(0,0) *= (1-damage);
    //     stiff(1,1) *= (1-damage);
    //   }
    //   if((B_matrix.transpose()*V_col)[1]>0.0) // Opening faces
    //   {
    //     stiff(0,0) *=  (1-damage);
    //     stiff(1,1) *= -1.0*damage_init/(damage_fini-damage_init);
    //
    //     if(abs(damage-1)<1E-8){
	  //        stiff(1,1)=0.0;
	  //       }
    //     // stiff(1,1) *= -1.0*damage_init/(damage_fini-damage_init);
    //   }
    // }

    MatrixXd Projection(2,2);  //Projects the separation onto the Axis
    Projection.setIdentity(); // Include at a later time

    Disp_stiff = B_matrix*(stiff*(Projection*B_matrix.transpose()));
};


void Matrices::Eval_LHS(float time, double* hht){

    LHS = Disp_stiff;

};

void Matrices::Eval_Rhs(const MatrixXd& Mass,const MatrixXd& U,const MatrixXd& DU,const MatrixXd& V,const MatrixXd& A,float time, float delta_t,
const VectorXd& shape_funcs, const int* sequence, double* hht, double* SVARS, int gp, const int lflag){

  MatrixXd U_col= U.reshaped<RowMajor>();
  MatrixXd DU_col= DU.reshaped<RowMajor>();
  MatrixXd V_col= V.reshaped<RowMajor>();
  MatrixXd A_col= A.reshaped<RowMajor>();
  VectorXd temp_SVARS = VectorXd::Zero(2*number_elements);

  // Update the state variables at GP
  // History variables - using semi group property of Exponential function
  //

  for(int i=0;i<2*number_elements;++i)
  {
    temp_SVARS(i) = *(2+SVARS+i+2*gp*number_elements);

  }

  for(int i=0;i<number_elements;++i)
  {
    temp_SVARS(seq(i*2,(i*2)+1))*= exp(-delta_t/cohesive_prony(i,1));
    temp_SVARS(seq(i*2,(i*2)+1))+= exp(-delta_t/2.0/cohesive_prony(i,1))*B_matrix.transpose()*DU_col;
  }

  // SVARS[(2+(i+gp*number_elements)*2)..(2+(i+gp*number_elements)*2)+1]
  // Update the state variables only when required
  if(lflag==1)
  {
    for(int i=0;i<2*number_elements;++i)
    {
      *(2+SVARS+i+2*gp*number_elements) = temp_SVARS(i) ;
    }
  }

  double alpha = *hht;
  double beta = *(hht+1);
  double gamma = *(hht+2);


  // Disp vector part
  MatrixXd stiff(2,2);
  stiff.setZero();

  stiff(0,0) = cohesive_stiff;
  stiff(1,1) = cohesive_stiff;

  if(damage>0)
  {
    stiff(0,0) *= (1-damage);
    stiff(1,1) *= (1-damage);
  }

  MatrixXd Projection(2,2);  //Projects the separation onto the nodal coordinate
  Projection.setIdentity();
  Rhs += -1.0*B_matrix*(stiff*(Projection*B_matrix.transpose()))*U_col - alpha*B_matrix*(stiff*(Projection*B_matrix.transpose()))*DU_col;
  cohesive_stiff /= (1.0-(cohesive_prony.colwise().sum())(0));
  for(int i=0;i<number_elements;++i)
  {
    Rhs += -1.0*(1-damage)*cohesive_prony(i,0)*cohesive_stiff*(B_matrix*temp_SVARS(seq(2*i,2*i+1)));
    // Rhs += -1.0*alpha*(1-damage)*cohesive_prony(i,0)*cohesive_stiff*(B_matrix*temp_SVARS(seq(2*i,2*i+1))        );
  }

};


void Matrices::evaluate_BMatrix(MatrixXd& B_matrix,const MatrixXd& shape_funcs,const int* sequence){

  B_matrix.resize(8,2);
  B_matrix(2*(*(sequence+0)),0) = -shape_funcs(0);
  B_matrix(2*(*(sequence+0))+1,0) = 0.0;
  B_matrix(2*(*(sequence+1)),0) = -shape_funcs(1);
  B_matrix(2*(*(sequence+1))+1,0) = 0.0;
  B_matrix(2*(*(sequence+2)),0) = shape_funcs(1);
  B_matrix(2*(*(sequence+2))+1,0) = 0.0;
  B_matrix(2*(*(sequence+3)),0) = shape_funcs(0);
  B_matrix(2*(*(sequence+3))+1,0) = 0.0;

  B_matrix(2*(*(sequence+0)),1) = 0.0;
  B_matrix(2*(*(sequence+0))+1,1) = -shape_funcs(0);
  B_matrix(2*(*(sequence+1)),1) = 0.0;
  B_matrix(2*(*(sequence+1))+1,1) = -shape_funcs(1);
  B_matrix(2*(*(sequence+2)),1) = 0.0;
  B_matrix(2*(*(sequence+2))+1,1) = shape_funcs(1);
  B_matrix(2*(*(sequence+3)),1) = 0.0;
  B_matrix(2*(*(sequence+3))+1,1) = shape_funcs(0);

};

// *********************************************************************************
//Access functions
MatrixXd Matrices::Get_LHS() const{
  return LHS;
};


MatrixXd Matrices::Get_Rhs() const{
  // // cout<<Rhs<<endl;
  return Rhs;
};
