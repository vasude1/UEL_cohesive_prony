/* Matrices CLASS TO DEFINE UELMAT ELMENT LHS AND RHS
Created by Va Ka
--------------------------------------------------------------------------------
Memory Management -

SVARS structure -----
d1(1),f11,f12,f13...f18(8),h11x,h11y,h12x,h12y,h13x,h13y(2*nel),....
d2(1),f21,f22,f23...f28(8),h21x,h21y,h22x,h22y,h31x,h31y(2*nel),....

SVARS 0 and 1 are the damage variables at Gauss points

24 terms at a gp in all

SVARS 2 to 9 contain the separation at previous time step at gp 1
SVARS 10 to 17 contain the force vector at gp 1 at time n-1
SVARS 18 to 19 contain the h-variable at gp 1
20 to 25 are backup

SVARS 26 to 33 contain the separation at previous time step at gp 2
SVARS 34 to 41 contain the force vector at gp 2 at time n-1
SVARS 42 to 43 contain the h-variable at gp 2
------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include<fstream>
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

  number_force = 8;
  number_elements = 6;
  number_h = 12;
  var_per_gp = 1+number_force+number_h;
  cohesive_prony = MatrixXd::Zero(number_elements,2);
  dt_gdt = 0.0;

  LHS.resize(ndofel,ndofel); // Gausspoint Stiffness contribution

  Disp_stiff.resize(ndofel,ndofel);

  Velo_stiff.resize(ndofel,ndofel);

  Rhs.resize(ndofel);

  nnodes = ndofel/2;


  cohesive_prony <<  0.9717,1E-7,
  0.01,  1e-06,
  0.00500388,  1e-05,
  0.00149131, 0.0001,
  0.00115114,  0.001,
  0.00121,   0.01;
  cohesive_prony(seq(0,5),0) *= 1.0/(1.0-cohesive_prony(0,0)-cohesive_prony(1,0)-cohesive_prony(2,0)-cohesive_prony(3,0)-cohesive_prony(4,0)-cohesive_prony(5,0));

  // std::cout << cohesive_prony << '\n';
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void Matrices::Set_Matrices(const int lflag,const ShapeFunc& Shape,const MatrixXd& Coords,const MatrixXd& U,const MatrixXd& DU,const MatrixXd& V,const MatrixXd& A,const MatrixXd& Mass ,double* SVARS,float delta_t,
                          float time, const int* sequence, int gp, double* hht)
{
  MatrixXd shape_funcs = Shape.Get_Shape();

  cohesive_stiff = 1E13;

  evaluate_BMatrix(B_matrix, shape_funcs, sequence);

  Eval_velocity_jump(B_matrix,V);

  define_damage_var(time,Coords);

  set_Zero();

  Eval_separation(shape_funcs,U,sequence);

  Eval_damage(lflag,separation,SVARS,gp);

  Eval_Rhs(Mass,U,DU,V,A,time,delta_t, shape_funcs, sequence, hht,SVARS,gp,lflag);

  Eval_stiff(separation, sequence, U, delta_t);

  Eval_LHS(time,hht);


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
    // if(x_centroid < 3E-3)  // Seed crack
    // {
    //   damage_init = 1E-10;
    //   damage_fini = 2E-3;
    // }
  if (time>1.0)   // time>1.0
  {

    if(x_centroid < 3E-3)  // Seed crack
    {
      damage_init = 1E-7;
      damage_fini = 2E-3;
    }
  }

};

void Matrices::Eval_velocity_jump(const MatrixXd& B_matrix, const MatrixXd& V){
  VectorXd V_col= V.reshaped<RowMajor>();
  VectorXd jump = B_matrix.transpose()*V_col;
  double jump_x = real(jump(0));
  double jump_y = real(jump(1));
  // velocity_jump = pow(pow(jump_x,2)+pow(jump_y,2),0.5);
  velocity_jump = jump_y;

};

void Matrices::Eval_separation(const VectorXd& shape_funcs,const MatrixXd& U, const int* sequence){

  //  MatrixXd left_separation;  // Separation on the Left
  // left_separation = U(*(sequence+3),seq(0,1)) - U(*(sequence+0),seq(0,1));
  //
  //  MatrixXd right_separation; //Separation on the right end
  // right_separation = U(*(sequence+2),seq(0,1)) - U(*(sequence+1),seq(0,1));
  //
  //  MatrixXd gp_separation;  // Separation at Gauss point = N1*u1+N2*u2
  // gp_separation = shape_funcs(0)*left_separation+shape_funcs(1)*right_separation;
  //
  // float separation_x = real(gp_separation(0,0));
  // float separation_y = real(gp_separation(0,1));
  // separation = pow(pow(separation_x,2)+pow(separation_y,2),0.5);

  VectorXd separation_vec = B_matrix.transpose()*U.reshaped<RowMajor>();
  // std::cout << U << '\n';
  // std::cout << "  " << '\n';
  separation = separation_vec.norm();
  if(separation_vec(1)<0)
  {
    separation = 0.0;
  }

};


void Matrices::Eval_damage(const int lflag, float separation,double* SVARS, int gp){

  // If fully damaged, return
  // std::cout << separation << '\n';
  damage=0.0;
  if(separation > damage_init)
  {
    // std::cout << separation << '\n';
    damage = damage_fini*(separation-damage_init)/(separation)/(damage_fini-damage_init);
  }

  if(damage > *(SVARS+gp*var_per_gp))
  {
    if(damage>=1)
    {
      damage=1.0;
    }
    if(lflag==1)
    {
        *(SVARS+gp*var_per_gp)=damage;
    }
    return;
  }

  if(damage < *(SVARS+gp*var_per_gp))
  {
    damage=*(SVARS+gp*var_per_gp);
    return;
  }
  return;

};


void Matrices::Eval_stiff(float separation,const int* sequence, const MatrixXd& U, float delta_t){

    // Elemental matrices  (1+dt_gdt)*
    // std::cout << separation << '\n';
    MatrixXd stiff = cohesive_stiff*MatrixXd::Identity(2,2);
    // dt_gdt = 0.0;
    // To determine if the faces are opening or closing
    // VectorXd V_col= V.reshaped<RowMajor>();
    // w = B_matrix.transpose()*U.reshaped<RowMajor>();
    Disp_stiff = (1-damage)*(1+dt_gdt)*B_matrix*(stiff*B_matrix.transpose());
    if((velocity_jump>0) && (separation>damage_init))
    {
      // std::cout << "Separation Exceeded" << '\n';
      Disp_stiff -= damage_init*damage_fini/(damage_fini-damage_init)/pow(separation,3)*B_matrix*stiff*(w*(B_matrix.transpose()*U.reshaped<RowMajor>()).transpose())*B_matrix.transpose();
    }
    // std::cout << 1+dt_gdt << '\n';
    if (damage>0) {
      std::ofstream ofile;
      ofile.open("/home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony/Res.txt", std::ios::app);
      ofile << "gdt = "<< '\n';
      ofile <<1+dt_gdt<< '\n';
      ofile << "Damage= "<<'\n';
      ofile << damage <<'\n';
      ofile << "Stiff= "<< '\n';
      ofile << Disp_stiff << '\n';
      ofile << " " <<'\n';
      ofile.close();
    }
};


void Matrices::Eval_LHS(float time, double* hht){

    LHS = (1+*hht)*Disp_stiff;
    // std::cout << LHS << '\n';

};

void Matrices::Eval_Rhs(const MatrixXd& Mass,const MatrixXd& U,const MatrixXd& DU,const MatrixXd& V,const MatrixXd& A,float time, float delta_t,
const VectorXd& shape_funcs, const int* sequence, double* hht, double* SVARS, int gp, const int lflag){
  // if (time>1.0 && separation > 1e-6)
  // {
  //   std::cout << separation << '\n';
  // }
  MatrixXd U_col= U.reshaped<RowMajor>();
  MatrixXd DU_col= DU.reshaped<RowMajor>();
  // MatrixXd V_col= V.reshaped<RowMajor>();
  // MatrixXd A_col= A.reshaped<RowMajor>();
  double alpha = *hht;
  double beta = *(hht+1);
  double gamma = *(hht+2);

  // To store the current force vector
  VectorXd f_n(number_force);

  // Compute Separations from displacements
  Vector2d dDelta_n = B_matrix.transpose()*DU_col;
  Vector2d Delta_n = B_matrix.transpose()*U_col;
  Vector2d Delta_nm1 = Delta_n - dDelta_n;

  VectorXd f_nm1(number_force);
  for(int i =0;i<number_force;++i)
  {
    f_nm1(i) = *(SVARS+gp*var_per_gp+i+1);
  }

  // Store Internal variables and old force
  VectorXd hn(number_h);
  for(int i=0;i<number_h;++i)
  {
    hn(i) = *(SVARS+gp*var_per_gp+number_force+i+1);
  }

  // Compute factors for force computation
  for(int i=0;i<number_elements;++i)
  {
    dt_gdt +=cohesive_prony(i,0)*exp(-delta_t/2.0/cohesive_prony(i,1));
  }

  // Compute the current force vector
  MatrixXd stiff = cohesive_stiff*MatrixXd::Identity(2,2);
  if(damage>0)
  {
    stiff = (1-damage)*stiff;
  }

  w = Delta_n + dt_gdt*dDelta_n;
  for(int i=0;i<number_elements;++i)
  {
    w += exp(-delta_t/cohesive_prony(i,1))*hn(seq(2*i,2*i+1));
    // std::cout << exp(-delta_t/cohesive_prony(i,1))*hn(seq(2*i,2*i+1)) << '\n';
  }
  f_n = B_matrix*stiff*w;
  // HHT = (1+alpha)*current force - alpha * old force
  Rhs = -1.0*(1+alpha)*f_n + alpha*f_nm1;

  if (damage>0) {
    std::ofstream ofile;
    ofile.open("/home/vasudevan/PhD/Code/UEL_cohesive/implicit/cohesive_prony/Res.txt", std::ios::app);
    ofile << "Delta_n= "<< '\n';
    ofile << Delta_n << '\n';
    ofile << "dDelta_n= "<< '\n';
    ofile << dDelta_n << '\n';
    ofile << "Rhs= "<< '\n';
    ofile << Rhs << '\n';
    ofile.close();
  }

  // lflag=1 update the internal variables
  if(lflag==1)
  {
    VectorXd temp_SVARS(2);

    for(int i=0;i<number_elements;++i)
    {
      temp_SVARS = cohesive_prony(i,0)*exp(-delta_t/2.0/cohesive_prony(i,1))*dDelta_n;
    *(SVARS+gp*var_per_gp+1+number_force+2*i) = (*(SVARS+gp*var_per_gp+number_force+2*i))*exp(-delta_t/cohesive_prony(i,1))+temp_SVARS(0) ;
    *(SVARS+gp*var_per_gp+1+number_force+2*i+1) = (*(SVARS+gp*var_per_gp+number_force+2*i+1))*exp(-delta_t/cohesive_prony(i,1))+temp_SVARS(1) ;
    }
    for(int i =0;i<number_force;++i)
    {
      *(SVARS+gp*var_per_gp+1+i) = f_n(i) ;
    }
  }

  if(lflag==5)
  {
    Rhs = 0.5*(f_n+f_nm1);
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
