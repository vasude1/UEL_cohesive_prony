
// Created by VaKa
// Shall be moved to Matrices class at a later time

#ifndef MASS_MAT_H
#define MASS_MAT_H

void shape(MatrixXd& shape_fns, float e, float n, int* order)
{
  shape_fns(2*(*(order+0)),0) = (1-e)*(1-n)/4.0;
  shape_fns(2*(*(order+1)),0) = (1+e)*(1-n)/4.0;
  shape_fns(2*(*(order+2)),0) = (1+e)*(1+n)/4.0;
  shape_fns(2*(*(order+3)),0) = (1-e)*(1+n)/4.0;

  shape_fns(2*(*(order+0))+1,1) = (1-e)*(1-n)/4.0;
  shape_fns(2*(*(order+1))+1,1) = (1+e)*(1-n)/4.0;
  shape_fns(2*(*(order+2))+1,1) = (1+e)*(1+n)/4.0;
  shape_fns(2*(*(order+3))+1,1) = (1-e)*(1+n)/4.0;
};



double computeJa(const MatrixXd verts,float e,float n){

  MatrixXd Ja(2,2);

  float x1=verts(0,0);
  float y1=verts(0,1);

  float x2=verts(1,0);
  float y2=verts(1,1);

  float x3=verts(2,0);
  float y3=verts(2,1);

  float x4=verts(3,0);
  float y4=verts(3,1);

  Ja(0,0)=(-x1+x2+x3-x4)/4+(x1-x2+x3-x4)*n/4;
  Ja(0,1)=(-y1+y2+y3-y4)/4+(y1-y2+y3-y4)*n/4;

  Ja(1,0)=((-x1-x2+x3+x4)/4)+((x1-x2+x3-x4)*e/4);
  Ja(1,1)=((-y1-y2+y3+y4)/4)+((y1-y2+y3-y4)*e/4);

  return Ja.determinant();

};

void compute_mass(MatrixXd& mass, const MatrixXd& coords, int* order)
{
  MatrixXd gauss_points(4,2);

  double usqt = 1.0/1.732;

  gauss_points(0,0) = -usqt;
  gauss_points(0,1) = -usqt;

  gauss_points(1,0) = usqt;
  gauss_points(1,1) = -usqt;

  gauss_points(2,0) = -usqt;
  gauss_points(2,1) = usqt;

  gauss_points(3,0) = usqt;
  gauss_points(3,1) = usqt;

  MatrixXd column(8,2);
  column.setZero();

  for(int gp=0; gp<4;++gp) // Loop over gauss points
  {
    float e = gauss_points(gp,0);
    float n = gauss_points(gp,1);
    shape(column,e,n, order);
    double Ja_det = computeJa(coords,e,n);
    mass += Ja_det*(column*column.transpose());
  }

  mass = 1044.0*mass;
};

 #endif
