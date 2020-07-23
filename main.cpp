/* TEST USER SUBROUTINE FOR UELMAT
Created by T.Tiirats // Date: 02/19
//--------------------------------*/

#include <stdio.h>
#include <iostream>
// #include <aba_for_c.h>
#include <Eigen/Dense> // Eigen class
#include <list>

#include "USUB_UtilityRoutines.h"

using namespace Eigen;
using namespace std;

#include "UELMAT_RTLM.cpp"

int main() {
	Utility Util(4,8,4);
	// 
	// double corrd [8] = {};
	// corrd[0] = 0.0;
	// corrd[1] = 0.0;
	// corrd[2] = 1.0;
	// corrd[3] = 0.0;
	// corrd[4] = 1.0;
	// corrd[5] = 1.0;
	// corrd[6] = 0.0;
	// corrd[7] = 1.0;
	//
	// double U [8] = {};
	// double DU [8] = {};
	//
	// U[0]=1.0;
	// U[1]=1.0;
	// U[2]=1.0;
	// U[3]=1.0;
	// U[4]=1.0;
	// U[5]=1.0;
	// U[6]=1.0;
	// U[7]=0.0;
	//
	//
	// 	DU[0]=1.0;
	// 	DU[1]=1.0;
	// 	DU[2]=1.0;
	// 	DU[3]=1.0;
	// 	DU[4]=1.0;
	// 	DU[5]=1.0;
	// 	DU[6]=1.0;
	// 	DU[7]=0.0;

	// for (int i=0; i<8;i++){
	// 	U[i]= 0.;
	// 	DU[i] = 0.;
	// }

	double AMATRX [64];
	double RHS[8];
	Util.Read_Input(corrd,U,DU);
	Util.Set_Output(AMATRX,RHS);
	UELMAT_RTLM(Util);
	// for (int i=0; i< 64;i++){
	// 	cout<<AMATRX[i]<<endl;
	// }
	// for (int i=0; i< 8;i++){
	// 	cout<<RHS[i]<<endl;
	// }
	// MatrixXd LHS1;
	// LHS1.resize(8,8);
	Map<MatrixXd> LHS1(AMATRX, 8,8);
	Map<MatrixXd> RHS1(RHS, 8,1);

	cout<<LHS1<<endl;
	cout<<RHS1<<endl;



	return 0;
}

// void UELMAT_RTLM(Utility& Util);
//
// extern "C"
// void FOR_NAME(uelmat)(double* RHS,double* AMATRX,double* SVARS,
// 		      double* ENERGY,int& NDOFEL,int& NRHS,int& NSVARS,
// 		      double* PROPS,int& NPROPS,double* COORDS,int& MCRD,
// 		      int& NNODE,double* U,double* DU,double* V,double* A,
// 		      int& JTYPE,double& TIME,double& DTIME,int& KSTEP,
// 		      int& KINC,int& JELEM,double* PARAMS,int& NDLOAD,
// 		      int& JDLTYP,double& ADLMAG,double* PREDEF,int& NPREDF,
// 		      int* LFLAGS,int& MLVARX,double& DDLMAG,int& MDLOAD,
// 		      double& PNEWDT,int& JPROPS,int& NJPROP,double& PERIOD,
// 		      double* MATERIALLIB){
//
//   cout << "START UELMAT ----- "<< endl;
//
//   // USER INPUT ________________________________________________________________
//   int integ_ord = 2; // Gauss quadrature order in tangential direction // to Properties
//
//
//   // READ VALUES FROM FORTRAN VARIABLES ________________________________________
//   Utility Util(NNODE,NDOFEL,integ_ord); // Utility routines to handle global variables
//   Util.Read_Input(COORDS,U,DU);            // Read Inputs
//   Util.Set_Output(AMATRX,RHS);          // Pass Pointers for output fields
//   list<int> lFlags( {LFLAGS[0],LFLAGS[1],LFLAGS[2],LFLAGS[3],LFLAGS[4]} );
//
//   // RUN ROUTINE (according to flags) __________________________________________
//   if (LFLAGS[2]==4){ // Define only mass matrix
//     Util.Set_Identity();
//     cout << " Return Mass Matrix (defined as Identity). " << endl;
//   }
//   else {
//     if (lFlags != list<int> {1,0,1,0,1} ){
// 	// Flag Check
// 	cout << "ERROR!!! Subroutine not suitable for this case. Check LFLAGS: " << endl;
// 	cout << " LFLAGS: "<<LFLAGS[0]<<" "<<LFLAGS[1]<<" "<<LFLAGS[2]<<" "
// 	     <<LFLAGS[3]<<" "<<LFLAGS[4]<<endl;
// 	throw exception();
//       };
//
//     // CALL UELMAT LIBRARY _____________________________________________________
//     UELMAT_RTLM(Util);
//   };
//
//   cout << "END UELMAT ----- "<< endl;
//   return;
// };
