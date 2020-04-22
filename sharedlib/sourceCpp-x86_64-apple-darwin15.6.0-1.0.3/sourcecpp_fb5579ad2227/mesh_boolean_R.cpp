#include <Rcpp.h>
#include <iostream>
#include <fstream>
//#include <cmath>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <Eigen/Core>
using namespace Rcpp;

  
  // [[Rcpp::export]]
List mesh_boolean( NumericMatrix VAi, NumericMatrix FAi, NumericMatrix VBi, NumericMatrix FBi,  std::string op_type) {
  

  int VA_nrow = VAi.nrow(), VA_ncol = VAi.ncol();
  int VB_nrow = VBi.nrow(), VB_ncol = VBi.ncol();
  int FA_nrow = FAi.nrow(), FA_ncol = FAi.ncol();
  int FB_nrow = FBi.nrow(), FB_ncol = FBi.ncol();
  
  // std::cout << VA_nrow << "\t" << VA_ncol << "\t" << VB_nrow << "\t" << VB_ncol << "\t" << FA_nrow << "\t" << FA_ncol << "\t" << FB_nrow << "\t" << FB_ncol << "\n" ; 
 

  Eigen::MatrixXd VA(VA_nrow,VA_ncol), VB(VB_nrow,VB_ncol),VC;
  Eigen::VectorXi J;
  Eigen::MatrixXi FA(FA_nrow,FA_ncol), FB(FB_nrow,FB_ncol),FC;

      for(int i=0; i < VA_nrow; i++){
        for(int j=0; j< VA_ncol; j++){
          VA(i,j) = VAi(i,j);
        }
      }

      for(int i=0; i < VB_nrow; i++){
        for(int j=0; j< VB_ncol; j++){
          VB(i,j) = VBi(i,j);
        }
      }

      for(int i=0; i < FA_nrow; i++){
        for(int j=0; j< FA_ncol; j++){
          FA(i,j) = FAi(i,j);
        }
      }

      for(int i=0; i < FB_nrow; i++){
        for(int j=0; j< FB_ncol; j++){
          FB(i,j) = FBi(i,j);
        }
      }
      
      if(op_type == "intersect"){
        igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_INTERSECT,VC,FC,J);
      }else if(op_type == "union"){
        igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_UNION,VC,FC,J);
      }else if(op_type == "xor"){
        igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_XOR,VC,FC,J);
      }else if(op_type == "minus"){
        igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_MINUS,VC,FC,J);
      }else if(op_type == "resolve"){
        igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_RESOLVE,VC,FC,J);
      }

      

      // std::cout << "a" << "\n";
      // std::cout << VC.rows() << "\t" << VC.cols() << "\t" << FC.rows() << "\t" << FC.cols() << "\n";
      
      NumericMatrix VCo(VC.rows(), VC.cols()), FCo(FC.rows(), FC.cols());

      for(int i=0; i < VC.rows(); i++){
        for(int j=0; j< VC.cols(); j++){
          VCo(i,j) = VC(i,j);
        }
      }

      for(int i=0; i < FC.rows(); i++){
        for(int j=0; j< FC.cols(); j++){
          FCo(i,j) = FC(i,j);
        }
      }


    return List::create(Named("VC") = VCo, Named("FC") = FCo);
    
    
  }

#include <Rcpp.h>
// mesh_boolean
List mesh_boolean(NumericMatrix VAi, NumericMatrix FAi, NumericMatrix VBi, NumericMatrix FBi, std::string op_type);
RcppExport SEXP sourceCpp_1_mesh_boolean(SEXP VAiSEXP, SEXP FAiSEXP, SEXP VBiSEXP, SEXP FBiSEXP, SEXP op_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type VAi(VAiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type FAi(FAiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type VBi(VBiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type FBi(FBiSEXP);
    Rcpp::traits::input_parameter< std::string >::type op_type(op_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(mesh_boolean(VAi, FAi, VBi, FBi, op_type));
    return rcpp_result_gen;
END_RCPP
}
