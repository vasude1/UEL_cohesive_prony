/* Matrices CLASS TO DEFINE UELMAT ELMENT LHS AND RHS
Created by VaKa
--------------------------------------------------------------------------------
Last edit:


------------------------------------------------------------------------------*/

// INCLUDES
#include <iostream>
#include <cmath>
#include <Eigen/Dense> // Eigen class

// OWN INCLUDES
#include "UELMAT_Order.h"

// Namespaces
using namespace Eigen;
using namespace std;

// Constructors ****************************************************************
Order::Order(const MatrixXd& coordinates){
  sortout(coordinates);
};

void Order::sortout(const MatrixXd& coordinates){ // Over written for now
  int order[4];
  // float x_center=0,y_center=0;
  //
  //
  // for(int i=0;i<4;++i){
  //     x_center += coordinates(i,0);
  //     y_center += coordinates(i,1);
  // }
  // x_center /= 4.0;
  // y_center /= 4.0;
  //
  // for(int i=0;i<4;++i){
  //   if (coordinates(i,0) < x_center && coordinates(i,1) < y_center){
  //     order[0] = i;
  //   }
  //   if (coordinates(i,0) > x_center && coordinates(i,1) < y_center){
  //     order[1] = i;
  //   }
  //   if (coordinates(i,0) > x_center && coordinates(i,1) > y_center){
  //     order[2] = i;
  //   }
  //   if (coordinates(i,0) < x_center && coordinates(i,1) > y_center){
  //     order[3] = i;
  //   }
  //
  // }
  order[0] = 2;
  order[1] = 3;
  order[2] = 0;
  order[3] = 1;
  for(int i=0;i<4;++i)
  {
      *(this->sequence+i) = *(order+i);
      // std::cout << *(order+i) << '\n';
  }

};
