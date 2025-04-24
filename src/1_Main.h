#if !defined(_Main_H)
#define _Main_H

#include <RcppArmadillo.h>
#include "2_Data.h"

class CMain {
	
 public:
		
  CMain(arma::vec Y0, int n_bins0, int no_grid_x0) ; // constructor
  ~CMain() ; // destructor
      
  void init() ; 
  arma::vec Print_Y() ; 
		
  double mBCV_obj_fn(double h0) ;	
	double IACT() ;
  arma::vec Print_LU_Bound_h() ; 
  void Read_LU_Bound_h(arma::vec LU_Bound_h_) ;
		    
 private:
		
  CData Data;
		
};

#endif  //_CMain_H
