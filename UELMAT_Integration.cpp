/* INTEGRATION CLASS TO DEFINE UELMAT ELMENT INTEGRATION
Created by T.Tiirats
----------------------------------------------------------------------
Last edit:


--------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Integration.h"

// Namespaces
using namespace Eigen;
using namespace std;




// CLASS INTEGRATION FUNCTIONS ///////////////////////////////////////


// Constructors ******************************************************
Integration::Integration(){
  points.resize(2,1);
  weights.resize(2,1);
  set_Gausspoints();
  set_weights();
};
//--------------------------------------------------------------------

Integration::set_Gausspoints(){
  float usqt = 1./sqrt(3.);
  // float usqt = 1.0;
  points << -usqt, usqt;
};
//--------------------------------------------------------------------

Integration::set_weights(){
  float one = 1.0;
  weights << one, one;
};

float Integration::Get_Gp(int number) const{
  return points(number);
};

float Integration::Get_weight(int number) const{
  return weights(number);
};
