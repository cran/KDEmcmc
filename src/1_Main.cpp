#include "1_Main.h"

RCPP_MODULE(cKDEmodule){
  
  using namespace R;
  using namespace Rcpp;
  
  class_<CMain>( "cKDE" )
    
    .constructor<arma::vec, int, int>()     // expose the default constructor
    
    .method("init", &CMain::init, "Initialize Y and empirical density")
    .method("Y", &CMain::Print_Y, "Print Y")	
    
    .method("mBCV_obj_fn", &CMain::mBCV_obj_fn, "mBCV_obj_fn")
    .method("IACT", &CMain::IACT, "IACT")		
		.property("LU_Bound_h", &CMain::Print_LU_Bound_h, &CMain::Read_LU_Bound_h ,"Lower and Upper Bound of h")
    
    ;
  
}       

////////////////////////////////////////
//
CMain::CMain(arma::vec Y0, int n_bins0, int no_grid_x0) {
	
  Data.Y = Y0 ; 
  Data.n_bins = n_bins0 ;  // no of bins for calculations for \sum_i \sum_j phi^{(4)}((Y_i-Y_j)/g) or delta_ij
  Data.no_grid_x = no_grid_x0 ;  // no of bins for numerical integration for zeta_{\hat{f}_h}(K_h)
  
  // Default values unless changed by Print_ statemetns
  Data.if_log_display = false ; 
  Data.init() ; 
  
}

CMain::~CMain(){ }

////////////////////////////////////////
void CMain::init() { 
  Data.init() ; 
}

arma::vec CMain::Print_Y() { 
  return Data.Y ; 
}

////////////////////////////////////////
double CMain::mBCV_obj_fn(double h0){
  return Data.mBCV_obj_fn(h0) ;
}

double CMain::IACT(){
  return Data.IAT_Y ;
}

arma::vec CMain::Print_LU_Bound_h(){
  return Data.LU_Bound_h ;
}
void CMain::Read_LU_Bound_h(arma::vec LU_Bound_h_){
  Data.LU_Bound_h = LU_Bound_h_ ;
}
