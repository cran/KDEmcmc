#include "2_Data.h"

CData::CData() {
}

//Destructor
CData::~CData(){
}

/////////////////////////////////

void CData::init(){
  
  // Y summary
   
  n_Y = Y.n_rows ; 
  sd_Y = stddev(Y) ; min_Y = Y.min() ; max_Y = Y.max() ; 
  
  // empirical density of |Y_i-Y_j|, i=1,...,n-1, j=i+1,...,n (n chooses 2) in discrete grid (histograms) 
  bin_value = arma::zeros<arma::vec>(n_bins+1) ; // Center of bin
  bin_count = arma::zeros<arma::vec>(n_bins+1) ; // 
  // ex.  |     |     |     |     | max-min |      when  min=-1  max=8  ->  max-min = 9
  //      0     2     4     6     8         10     n_bins = 5 ;  bin_width = (4-(-4)) / (5-1) = 2
  
  // Replace approximaiton technique in SEXP bw_den(SEXP nbin, SEXP sx) of bandwdiths.c 
  // of stats in R base package
  /* This would be impracticable for long vectors.  Better to bin x first */
  
  bin_width = (max_Y-min_Y) * 1.01 / n_bins ;
  for (int i=1; i<n_bins+1; i++){
    bin_value(i) = bin_value(i-1) + bin_width ; 
  }
  
  for (int i = 0; i < (n_Y-1); i++) {
    int ii = floor(Y(i) / bin_width);
    // int ii = (int)(Y[i] / dd);	// Confirmed that this original code is wrong
    for (int j = (i+1); j < (n_Y); j++) {
      int jj = floor(Y(j) / bin_width);
      bin_count(abs(ii - jj)) = bin_count(abs(ii - jj)) + 1 ;
    }
  }
  
  //
  
  grid_x = arma::vec(no_grid_x) ;
  grid_x(0) = mean(Y) - 6.0 * stddev(Y) ; // minimum value
  meshsize_grid_x = 12.0 * stddev(Y) / (no_grid_x-1)  ; 
  for (int i_temp=1; i_temp<no_grid_x; i_temp++){
    grid_x(i_temp) = grid_x(i_temp-1) + meshsize_grid_x ;
  }
  
  // IACT
  IAT_Y = IACT_fn(Y) ; 
  
  // Lower and Upper Bound of h
  arma::vec sorted_Y = sort( Y ) ; // sort in "ascend" direction 
  double loc_25 = floor(n_Y * 0.25 ) ; double loc_75 = floor(n_Y * 0.75 ) ;
  double IQR = sorted_Y(loc_75) - sorted_Y(loc_25) ; 
  // std::cout << "25% of n=" << loc_25 << ", 75% of n=" << loc_75 << ", IQR=" << IQR << std::endl ; 
  
  LU_Bound_h = arma::vec(2) ; 
  LU_Bound_h(1) = IQR * pow(n_Y,(-1.0/5.0)) * pow(IAT_Y,(1.0/5.0)) ; 
  // std::cout << "n_Y=" << n_Y << ", pow(n_Y,(-1.0/5.0))=" << pow(n_Y,(-1.0/5.0)) << ", pow(IAT_Y,(1/5))=" << pow(IAT_Y,(1/5)) << std::endl ; 
  LU_Bound_h(0) = 0.1 * LU_Bound_h(1) ; 
  
  // NOTE
  // bw.SJ of R stat package uses   hmax = min( 1.144 * sd(x) * n^(-1/5) , 1.144 * IQR(x)/1.349 * n^(-1/5) )
  // bw.bcv of R stat package uses  hmax <- 1.144 * sd(x) * n^(-1/5)
  
} // void CData::init()

/////////////////////////////////

double CData::mBCV_obj_fn(double h0){
  
  double sum_divided_by_n = 0.0 ; 
  double Delta_ij, Delta_ij_2, Delta_ij_4 ; 
  
  for (int i=0; i<n_bins+1; i++){
    
    Delta_ij = bin_value(i) / h0; 
    Delta_ij_2 = Delta_ij * Delta_ij ; 
    Delta_ij_4 = Delta_ij_2 * Delta_ij_2 ; 
    
    if (Delta_ij_2 >= 1000){
      break;
    }
    // Then, exp(-1/4 Delta_ij^2) = 0 // Avoid slow and possibly error-producing underflows by cutting off at plus/minus sqrt(DELMAX) std deviations // Formulae (6.67) and (6.69) of Scott (1992), the latter corrected.
    
    sum_divided_by_n = sum_divided_by_n + ( 1.0 / n_Y ) * ( Delta_ij_4 - 12 * Delta_ij_2 + 12 ) * exp(-1.0 / 4.0 * Delta_ij_2) * bin_count(i) ;
  }
  
  double R_K = 1.0 / (2.0 * sqrt(M_PI)) ;
  double mbcv_obj_value =  1.0 / ( n_Y * h0) * R_K * zeta_K_h(h0) + 1.0 / (64.0 * n_Y * h0 * sqrt(M_PI)) * sum_divided_by_n ; // MODIFIED for large n_X
  
  return mbcv_obj_value ; 
  
}

/////////////////////////////////

double CData::IACT_fn(arma::vec y_vec){

  int n = y_vec.n_rows ;
  double y_bar = 1.0 / n * sum(y_vec) ;

  double gamma0 = 0.0 ; // variance estimator
  for (int i=0; i<n; i++){
    gamma0 = gamma0 + 1.0 / n * (y_vec(i)-y_bar)*(y_vec(i)-y_bar)  ;  // biased variance estimator with 1/n
  }

  int maxlag = floor(n/2) ;
  if ( maxlag <= 3 ) {
    Rcpp::stop("Not enough data, floor(n/2) <= 3 \n") ;  // break ;
  } ;

  ////////
  double gamma1 = 0.0 ; int lag = 1 ;
  for (int i=0; i<(n-lag); i++){
    gamma1 = gamma1 + 1.0 / n * ( y_vec(i)-y_bar )*( y_vec(i+lag)-y_bar ) ;
  }
  double old_Gamma = gamma0 + gamma1 ;

  double IAT = 1 + 2.0 * gamma1 / gamma0 * ( 1 - 1.0/n ) ; // Ver 1.6.2


  ////////
  int M = 1 ;
  double gamma_2M = 0.0 ; int lag_2M = 2 * M ;  // Ver 1.6.2
  for (int i=0; i<(n-lag_2M); i++){  // Ver 1.6.2
    gamma_2M = gamma_2M + 1.0 / n * ( y_vec(i)-y_bar )*( y_vec(i+lag_2M)-y_bar ) ;  // Ver 1.6.2
  }
  double gamma_2M1 = 0.0 ; int lag_2M1 = 2 * M + 1 ;  // Ver 1.6.2
  for (int i=0; i<(n-lag_2M1); i++){  // Ver 1.6.2
    gamma_2M1 = gamma_2M1 + 1.0 / n * ( y_vec(i)-y_bar )*( y_vec(i+lag_2M1)-y_bar ) ;  // Ver 1.6.2
  }
  double new_Gamma = gamma_2M + gamma_2M1 ;  // Ver 1.6.2

  ////////
  while ( (new_Gamma > 0) && (new_Gamma < old_Gamma) ) {

    IAT = IAT + 2.0 * ( (1.0 - 1.0*lag_2M/n)*gamma_2M + (1.0 - 1.0*lag_2M1/n)*gamma_2M1 ) / gamma0 ;  // Ver 1.6.2
    old_Gamma = new_Gamma ;  // Ver 1.6.2

    M = M + 1 ;
    if (2 * M + 1 > maxlag) {
      Rcpp::stop("Not enough data, maxlag=", maxlag, "\n") ;
      // break ;
    }

    gamma_2M = 0.0 ; lag_2M = 2 * M ;  // Ver 1.6.2
    for (int i=0; i<(n-lag_2M); i++){  // Ver 1.6.2
      gamma_2M = gamma_2M + 1.0 / n * ( y_vec(i)-y_bar )*( y_vec(i+lag_2M)-y_bar ) ;  // Ver 1.6.2
    }
    gamma_2M1 = 0.0 ; lag_2M1 = 2 * M + 1 ;  // Ver 1.6.2
    for (int i=0; i<(n-lag_2M1); i++){  // Ver 1.6.2
      gamma_2M1 = gamma_2M1 + 1.0 / n * ( y_vec(i)-y_bar )*( y_vec(i+lag_2M1)-y_bar ) ;  // Ver 1.6.2
    }
    new_Gamma = gamma_2M + gamma_2M1 ;

  } // while ( (new_Gamma > 0) && (new_Gamma < old_Gamma) )

  return IAT ;

} // double CData::IACT_fn(arma::vec y_vec)

/////////////////////////////////

double CData::zeta_K_h(double h0){

  arma::vec K_h(n_Y) ;

  double result_sum = 0.0 ;
  for (int i_point=0; i_point<no_grid_x; i_point++){
    for (int i_n=0; i_n<n_Y; i_n++){
      K_h(i_n) = 1.0/h0 * R::dnorm( (grid_x(i_point)-Y(i_n))/h0, 0, 1, 0) ; // log = F
    }
    double point_eval_f_hat = 1.0 / n_Y * sum(K_h) ;
    if ( point_eval_f_hat > 0 ) {
      // If all K_h are zeros, IACT_fn(K_h) = NaN
      // double tau_K_hx = IACT_fn(K_h) ;
      double tau_K_hx = IACT_for_tau_Khx_fn(K_h) ; 	// modified Version: 1.6.1, Date: 2017-10-28
      if ( std::isnan(tau_K_hx)==FALSE ){
        result_sum = result_sum + tau_K_hx * point_eval_f_hat * meshsize_grid_x ;
      }
    }
  }

  return result_sum ;

}

/////////////////////////////////

double CData::IACT_for_tau_Khx_fn(arma::vec K_h_vec){

  int n = K_h_vec.n_rows ;
  double K_h_bar = 1.0 / n * sum(K_h_vec) ;

  double gamma0 = 0.0 ; // variance estimator
  for (int i=0; i<n; i++){
    gamma0 = gamma0 + 1.0 / n * (K_h_vec(i)-K_h_bar)*(K_h_vec(i)-K_h_bar)  ;  // biased variance estimator with 1/n
  }

  int maxlag = floor(n/2) ;
  if ( maxlag <= 3 ) {
    Rcpp::stop("Not enough data, floor(n/2) <= 3 \n") ;  // break ;
  } ;

  ////////
  double gamma1 = 0.0 ; int lag = 1 ;
  for (int i=0; i<(n-lag); i++){
    gamma1 = gamma1 + 1.0 / n * ( K_h_vec(i)-K_h_bar )*( K_h_vec(i+lag)-K_h_bar ) ;
  }
  double old_Gamma = gamma0 + gamma1 ;

  double IACT_for_tau = 1 + 2.0 * gamma1 / gamma0 * ( 1 - 1.0/n ) ; // Ver 1.6.2

  ////////
  int M = 1 ;
  double gamma_2M = 0.0 ; int lag_2M = 2 * M ;  // Ver 1.6.2
  for (int i=0; i<(n-lag_2M); i++){  // Ver 1.6.2
    gamma_2M = gamma_2M + 1.0 / n * ( K_h_vec(i)-K_h_bar )*( K_h_vec(i+lag_2M)-K_h_bar ) ;   // Ver 1.6.2
  }
  double gamma_2M1 = 0.0 ; int lag_2M1 = 2 * M + 1 ;
  for (int i=0; i<(n-lag_2M1); i++){  // Ver 1.6.2
    gamma_2M1 = gamma_2M1 + 1.0 / n * ( K_h_vec(i)-K_h_bar )*( K_h_vec(i+lag_2M1)-K_h_bar ) ;   // Ver 1.6.2
  }
  double new_Gamma = gamma_2M + gamma_2M1 ;   // Ver 1.6.2

  ////////
  while ( (new_Gamma > 0) && (new_Gamma < old_Gamma) ) {

    IACT_for_tau = IACT_for_tau + 2.0 * ( (1.0 - 1.0*lag_2M/n)*gamma_2M + (1.0 - 1.0*lag_2M1/n)*gamma_2M1 ) / gamma0 ;  // Ver 1.6.2
    old_Gamma = new_Gamma ;  // Ver 1.6.2

    M = M + 1 ;
    if (2 * M + 1 > maxlag) {
      Rcpp::stop("Not enough data for IACT_for_tau_Khx, maxlag=", maxlag, "\n") ;
      // break ;
    }

    gamma_2M = 0.0 ; lag = 2 * M ;  // Ver 1.6.2
    for (int i=0; i<(n-lag_2M); i++){  // Ver 1.6.2
      gamma_2M = gamma_2M + 1.0 / n * ( K_h_vec(i)-K_h_bar )*( K_h_vec(i+lag_2M)-K_h_bar ) ;   // Ver 1.6.2
    }
    gamma_2M1 = 0.0 ; lag = 2 * M + 1 ;  // Ver 1.6.2
    for (int i=0; i<(n-lag_2M1); i++){  // Ver 1.6.2
      gamma_2M1 = gamma_2M1 + 1.0 / n * ( K_h_vec(i)-K_h_bar )*( K_h_vec(i+lag_2M1)-K_h_bar ) ;   // Ver 1.6.2
    }
    new_Gamma = gamma_2M + gamma_2M1 ;   // Ver 1.6.2

  } // while ( (new_Gamma > 0) && (new_Gamma < old_Gamma) )

  return IACT_for_tau ;

} // CData::IACT_for_tau_Khx_fn(arma::vec K_h_vec)
