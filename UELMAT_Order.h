/* HEADER FOR Matrices CLASS
Created by T.Tiirats
------------------------------------------------------------------------------*/


#ifndef  UELMAT_ORDER_H
#define  UELMAT_ORDER_H

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>// Eigen class

using namespace Eigen;
using namespace std;

class Order
{
public:
  int sequence[4];


//----------------------------------------------------------------------------

  // PUBLIC FUNCTIONS //////////////////////////////////////////////////////////
 public:
  // Signatures
  Order(const MatrixXd&);
  void sortout(const MatrixXd&);

};
#endif
