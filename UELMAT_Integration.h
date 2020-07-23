/* HEADER FOR INTEGRATION CLASS
Created by T.Tiirats
-----------------------------*/

#ifndef  UELMAT_Integration_H
#define  UELMAT_Integration_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

using namespace Eigen;
using namespace std;

class Integration
{
  // INTERNAL PARAMETERS /////////////////////////////////////////////
  int number_points = 2;  // Integration order
  MatrixXd points;
  MatrixXd weights;
 public:
  // Construct Integration schemes
  Integration();        // X-Directional Gauss
  set_Gausspoints();
  set_weights();
  float Get_Gp(int number) const;
  float Get_weight(int number) const;
};
#endif
