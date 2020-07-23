/* UELMAT SUBROUTINE FOR RTLM
Created by VaKa
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Integration.h"
#include "UELMAT_ShapeFunctions.h"
#include "UELMAT_Material.h"
#include "UELMAT_Assembly.h"
#include "USUB_UtilityRoutines.h"
#include "UELMAT_Order.h"


// Namespaces
using namespace Eigen;
using namespace std;

#include "UELMAT_mass_mat.h" //For mass matrix

void UELMAT_RTLM(Utility& Util, int lflag)
{

  MatrixXd Coords = Util.Coords;    // Nodal Coordinates
  MatrixXd U = Util.u;    // Nodal Displacements
  MatrixXd DU = Util.du;    // Nodal Displacement increments
  MatrixXd V = Util.v;    // Nodal Velocities
  MatrixXd A = Util.a;    // Nodal Accelerations
  float dt = Util.dt;         // Time_increment
  float TIME = Util.time;     // Current time
  double* _SVAR = Util.SVAR;  // POinter to state variables
  double* integ_para = Util.parameters; // Time integration parameters


  MatrixXd K;       //Stiffness Matrix
  VectorXd R;       //Residual Vector

  int problem_dimension = 1;  // One dimension cohesive element - zero thickness Remove later
  int ndofel = 8;             // 4 nodes times 2 DOF per node

  // MAss matrix conputation
  MatrixXd Mass(8,8);
  Mass.setZero();


  Order order(Coords);        // Get the sequence of nodal coordinates - connectivity

  Integration Integ;    // Integration rule for X direction

  ShapeFunc Shape(2,problem_dimension);    // Shapefunction manager

  Matrices Matrix(ndofel);    // Class computing stiffness and residuals

  Assembly Assmbl(ndofel);    // Gauss point assemble manager

  // compute_mass(Mass, Coords, order.sequence); // Caution - Valid only for zero thickness cohesive elements


  // Loop over the gauss points - Numerical integration

  for (int it_gp=0; it_gp < 2; ++it_gp){   // Two gauss points only

    Shape.Set_Approx(Integ.Get_Gp(it_gp), Coords, U, order.sequence);        // Compute shape functions at Gauss points

    Matrix.Set_Matrices(lflag,Shape,Coords,U,DU,V,A,Mass,_SVAR, dt, TIME, order.sequence ,it_gp, integ_para); // Compute the K and R

    float jac_det = Shape.Get_Jacobian_det(); // In this case, this is l/2 - 1D element

    Assmbl.Eval_Assembly(Matrix.Get_LHS(), Matrix.Get_Rhs(),
                          Integ.Get_weight(it_gp), jac_det);  // Assemble the matrices

  };


  switch (lflag) {
    case 1: //LFLAGS(3) = 1 Jacobian for AMATX and Residual for RHS
    {
      K = Assmbl.Get_Tangent();
      R = Assmbl.Get_Res();
      break;
    }

    case 2: //LFLAGS(3) = 2 Jacobian for AMATX
    {
      K = Assmbl.Get_Tangent();
      R.setZero();
      break;
    }

    case 4: //LFLAGS(3) = 4 Mass for AMATX
    {
      K = Mass;
      R.setZero();
      break;
    }

    case 5: //LFLAGS(3) = 5 Residual for RHS
    {
      K.setZero();
      R = Assmbl.Get_Res();
      break;
    }

    case 6: //LFLAGS(3) = 6 Mass for AMATX and Residual for RHS
    {
      K = Mass;
      R = Assmbl.Get_Res();
      break;
    }
  }

  Util.Write_Matrix(K,K.rows(),K.cols() );  // AMATRX
  Util.Write_Vector(R,R.size());                  // RHS

}
