#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/volume.h>
#include <Eigen/Core>
using namespace Rcpp;

  
  // [[Rcpp::export]]
double mesh_volume(NumericMatrix Vi, NumericMatrix Fi) {
  

  int V_nrow = Vi.nrow(), F_nrow = Fi.nrow();
  int V_ncol = Vi.ncol(), F_ncol = Fi.ncol();
  
  // std::cout << VA_nrow << "\t" << VA_ncol << "\t" << VB_nrow << "\t" << VB_ncol << "\t" << FA_nrow << "\t" << FA_ncol << "\t" << FB_nrow << "\t" << FB_ncol << "\n" ; 
 

  Eigen::MatrixXd V(V_nrow,V_ncol);
  Eigen::MatrixXi F(F_nrow,F_ncol);

  for(int i=0; i < V_nrow; i++){
    for(int j=0; j< V_ncol; j++){
      V(i,j) = Vi(i,j);
    }
  }

  
  for(int i=0; i < F_nrow; i++){
    for(int j=0; j< F_ncol; j++){
      F(i,j) = Fi(i,j);
    }
  }

  Eigen::MatrixXd V2(V.rows() + 1, V.cols());
  V2.topRows(V.rows()) = V;
  V2.bottomRows(1).setZero();
  Eigen::MatrixXi T(F.rows(), 4);
  T.leftCols(3) = F;
  T.rightCols(1).setConstant(V.rows());
  Eigen::VectorXd vol;
  igl::volume(V2, T, vol);

  return std::abs( vol.sum());
  }