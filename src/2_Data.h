#if !defined(_CData_H)
#define _CData_H

#include <RcppArmadillo.h>

class CData {

 public:
  
  CData(); // constructor
  ~CData(); // destructor
  
  void init() ; 
  
  //
  bool if_log_display ;  
  arma::vec Y ; 
	
  int n_bins ;
  arma::vec bin_value, bin_count ;
  double bin_width ; 

  double IACT_fn(arma::vec Y0) ;  
  double IACT_for_tau_Khx_fn(arma::vec K_h_vec)  ; 

  arma::vec grid_x ; 
  double meshsize_grid_x ; 
  int no_grid_x ; 
    
  double zeta_K_h(double h0) ; 
  
  arma::vec LU_Bound_h ; 
 
  int n_Y ; 
  double sd_Y, min_Y, max_Y ; 

  double mBCV_obj_fn(double h0) ;	
  double IAT_Y ;
	  
};

#endif  //_CData_H
