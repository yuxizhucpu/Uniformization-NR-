#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat matrixExp(arma::mat M) {
  // get function from environment
  Function expm("expm");
  // compute matrix exponential
  NumericMatrix res = expm(wrap(M),"Pade");
  // convert back to arma::mat and return
  return as<arma::mat>(res);
}

namespace msm {

List matrixPower(arma::mat mat, int power) {
  
  List to_return(power);
  arma::mat mat1;
  to_return[0] = mat1.eye(3,3);
  to_return[1] = mat;
  arma::mat result = mat;
  for(int i = 2; i < power; i++) {
    result = result * mat;
    to_return[i]=result;
  }
  return to_return;
}


List matrixPowerfd(arma::mat mat1,arma::mat qd, int power) {
  arma::mat pd = qd;
  
  List to_return(power); 
  to_return[0] = qd;
  List powerstore = msm::matrixPower(mat1,power); 
  
  for(int i = 1; i <  power; i++) {
    
    arma::mat stepi = powerstore[i-1];
    
    pd = stepi * qd+pd * mat1;
    to_return[i]=pd;
  }
  
  return to_return;
  
}

long double factorial(const int n)
{
  long f = 1;
  for (int i=1; i<=n; ++i)
    f *= i;
  return f;
}

long double logfactorial(const int n)
{
  long double f = 0;
  for (int i=1; i<=n; ++i)
    f += log(i);
  return f;
}

arma::mat padeApproximation(const arma::mat& A, int degree) {
  int n = A.n_rows;
  arma::mat identity = arma::eye<arma::mat>(n, n);
  arma::mat scaledA = A / pow(2, degree);
  arma::mat scaledA2 = scaledA * scaledA;
  
  arma::mat U = (identity + scaledA);
  arma::mat V = (identity - scaledA);
  arma::mat P = arma::zeros<arma::mat>(n, n);
  arma::mat Q = arma::zeros<arma::mat>(n, n);
  
  P.row(0) = V.row(0);
  Q.row(0) = U.row(0);
  
  for (int k = 1; k <= degree; ++k) {
    P = (scaledA2 * P);
    Q = (scaledA2 * Q);
    P += V.row(k) * Q(0, 0);
  }
  
  arma::mat result = U.col(0) * P;
  result = result * inv(Q);
  
  for (int i = 1; i <= degree; ++i) {
    result *= (identity + 2 * P / Q(0, 0));
  }
  
  return result;
}


}

// [[Rcpp::export]]
double likelihood(arma::mat& mat1) {
  std::map<int, double> timeMap;
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat jump_p=p1;
  
    // Create an identity matrix
    arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
    arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
  
    double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
      std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);

    if (timeMap.count(msmform(16)) == 0) { 
      
      for (int j=0; j < 45; j++) {
        
        double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * std::log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
        jump_p=jump_p * p1;
        
        // Convert matrix to vector by row
        jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
        kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
  
      }
      to_return_final += -log(kthpvector);
      timeMap[msmform(16)]= -log(kthpvector);
      
    } else {
      to_return_final += timeMap[msmform(16)];
    }
    
   
  }
  
  return to_return_final;
  }

// [[Rcpp::export]]
double likelihoodno_map(arma::mat& mat1) {

  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat jump_p=p1;
    
    // Create an identity matrix
    arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
    arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
    
    double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
      std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);
    
      for (int j=0; j < 45; j++) {
        
        double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * std::log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
        jump_p=jump_p * p1;
        
        // Convert matrix to vector by row
        jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
        kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
       
    }
      to_return_final += -log(kthpvector);
    
  }
  
  return to_return_final;
}

/*use map too slow as the data is small*/
// [[Rcpp::export]]
double likelihood_ward(arma::mat& mat1) {
  
  std::map<int, double> timeMap;
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    if (timeMap.count(msmform(11)) == 0) { 
      // Convert vector to matrix by row
      int numRows = 3;  // Number of rows in the resulting matrix
      int numCols = 3 ;  // Number of columns in the resulting matrix
      arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();

      arma::mat trans_p = arma::expmat( p1* msmform(9));
      
      arma::vec trans_p_vector = arma::vectorise(trans_p, 1).t();
      
      to_return_final += -log(trans_p_vector(msmform(10)-1));
      
      timeMap[msmform(11)]=-log(trans_p_vector(msmform(10)-1));
      
    } else {
      to_return_final += timeMap[msmform(11)];
    }
  
  }
  
  return to_return_final;
}

// [[Rcpp::export]]
double likelihood_ward_no_map(arma::mat& mat1) {
  
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat trans_p = arma::expmat( p1* msmform(9));
      
    arma::vec trans_p_vector = arma::vectorise(trans_p, 1).t();
      
    to_return_final += -log(trans_p_vector(msmform(10)-1));
  }
  return to_return_final;
}

//use expm 
// [[Rcpp::export]]
double likelihood_ward2(arma::mat& mat1) {
  
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat trans_p =  matrixExp( p1*msmform(9));
    
    arma::vec trans_p_vector = arma::vectorise(trans_p, 1).t();
    
    to_return_final += -log(trans_p_vector(msmform(10)-1));
  }
  
  return to_return_final;
}

// [[Rcpp::export]]
double likelihood_decompno_map(arma::mat& mat1) {

  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();

      // Compute the eigenvalue decomposition
      arma::cx_vec eigval;
      arma::cx_mat eigvec;
      arma::eig_gen(eigval, eigvec, p1);
      
      // Compute the exponential of the eigenvalues
      arma::cx_vec expEigval = arma::exp(eigval*msmform(9));
      
      // Compute the matrix exponential
      arma::cx_mat expA = eigvec * arma::diagmat(expEigval) * arma::inv(eigvec);
      
      // Convert the complex matrix to a real-valued vector
      arma::vec trans_p_vector = arma::vectorise(arma::real(expA), 1).t();
      
      
      
      to_return_final += -std::log(trans_p_vector(msmform(10)-1));
      
     
      

    
  }
  
  return to_return_final;
}

// [[Rcpp::export]]
double likelihood_decomp(arma::mat& mat1) {
  std::map<int, double> timeMap;
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    if (timeMap.count(msmform(11)) == 0) { 
      
      // Compute the eigenvalue decomposition
      arma::cx_vec eigval;
      arma::cx_mat eigvec;
      arma::eig_gen(eigval, eigvec, p1);
      
      // Compute the exponential of the eigenvalues
      arma::cx_vec expEigval = arma::exp(eigval*msmform(9));
      
      // Compute the matrix exponential
      arma::cx_mat expA = eigvec * arma::diagmat(expEigval) * arma::inv(eigvec);
      
      // Convert the complex matrix to a real-valued vector
      arma::vec trans_p_vector = arma::vectorise(arma::real(expA), 1).t();
      
      
      
      to_return_final += -std::log(trans_p_vector(msmform(10)-1));
      
      timeMap[msmform(11)]=-log(trans_p_vector(msmform(10)-1));
      
    } else {
      
      to_return_final += timeMap[msmform(11)];
    }
    
  
  }
  
  return to_return_final;
}

// [[Rcpp::export]]
double likelihood_taylor(arma::mat& mat1) {
  
  int nrow_Mat = mat1.n_rows;
  double  to_return_final=0.0;
  
  // Create an identity matrix
  arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
  arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    const arma::vec& msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    arma::mat jump_p=p1;
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
    
    double kthpvector = identity_mat1_vector(msmform(10)-1)+ jump_p_vector(msmform(10)-1)*msmform(9) ;
    
    for (int j=0; j<20; j++) {
      
      double timeppvector =std::exp((j+2) * log(msmform(9)) - (msm::logfactorial(j+2)));
      jump_p=jump_p * p1;
      // Convert matrix to vector by row
      jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
      
      kthpvector += jump_p_vector(msmform(10)-1)*timeppvector;
      
    }

    to_return_final += -std::log(kthpvector);
  }
  
  return to_return_final;
}

// [[Rcpp::export]]
std::vector<double> likelihood_part(arma::mat mat1) {
  
  std::map<int, std::vector<double>> timeMap;
  
  int nrow_Mat = mat1.n_rows;
  std::vector<double>  to_return_final(9);
  
  double kthdvector_all  =0;
  double kthdvector1_all =0;
  double kthdvector2_all =0;
  
  double kthd2vector_all =0;
  double kthd2vector1_all=0;
  double kthd2vector2_all=0;
  double kthd2vector3_all=0;
  double kthd2vector4_all=0;
  double kthd2vector5_all=0;
  
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    arma::vec msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat jump_p=p1;
    
    arma::mat first_de1 = {{-1/msmform(12), 1/msmform(12), 0},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de2 = {{-1/msmform(12), 0, 1/msmform(12)},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de3 = {{0, 0, 0},
    {0, -1/msmform(12), 1/msmform(12)},
    {0, 0, 0}};
    
    
    
    arma::mat kthd1=first_de1;
    arma::mat kthd2=first_de2;
    arma::mat kthd3=first_de3;
    
    arma::mat kthd22 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_1 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_2 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_3 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_4 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_5 = arma::zeros<arma::mat>(3, 3);
    
    
    arma::mat jump_initial = arma::eye<arma::mat>(3, 3);
    
    
    // Create an identity matrix
    arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
    arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
    
    arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
    arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
    arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
    
    
    double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
      std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);
    
    
    
    double kthdvector = kthd1_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector1 = kthd2_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector2 = kthd3_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    
    double kthd2vector  = 0;
    double kthd2vector1 = 0;
    double kthd2vector2 = 0;
    double kthd2vector3 = 0;
    double kthd2vector4 = 0;
    double kthd2vector5 = 0;
    
    if (timeMap.count(msmform(16)) == 0) {
      
      for (int j=0; j<20; j++) {
        //double timeppvector =std::exp(-msmform(12) * msmform(13)) * (std::pow(msmform(12) * msmform(13), j + 2)) / msm::factorial(j + 2);
        // numerical flow;
        double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
        
        kthd22 = kthd1 * first_de1 + kthd22 * p1 + kthd1 * first_de1 ;
        kthd22_1 = kthd2 * first_de1 + kthd22_1 * p1 + kthd1 * first_de2;
        kthd22_2 = kthd3 * first_de1 + kthd22_2 * p1 + kthd1 * first_de3;
        
        kthd22_3 = kthd2 * first_de2 + kthd22_3 * p1 + kthd2 * first_de2;
        kthd22_4 = kthd3 * first_de2 + kthd22_4 * p1 + kthd2 * first_de3;
        
        kthd22_5 = kthd3 * first_de3 + kthd22_5 * p1 + kthd3 * first_de3;
        
        arma::vec kthd22_vector = arma::vectorise(kthd22, 1).t();
        arma::vec kthd22_1_vector = arma::vectorise(kthd22_1,1).t();
        arma::vec kthd22_2_vector = arma::vectorise(kthd22_2, 1).t();
        arma::vec kthd22_3_vector = arma::vectorise(kthd22_3, 1).t();
        arma::vec kthd22_4_vector = arma::vectorise(kthd22_4, 1).t();
        arma::vec kthd22_5_vector = arma::vectorise(kthd22_5, 1).t();
        
        kthd2vector += kthd22_vector(msmform(14)-1) * timeppvector;
        kthd2vector1 += kthd22_1_vector(msmform(14)-1) * timeppvector;
        kthd2vector2 += kthd22_2_vector(msmform(14)-1) * timeppvector;
        kthd2vector3 += kthd22_3_vector(msmform(14)-1) * timeppvector;
        kthd2vector4 += kthd22_4_vector(msmform(14)-1) * timeppvector;
        kthd2vector5 += kthd22_5_vector(msmform(14)-1) * timeppvector;
        
        kthd1 = jump_p * first_de1 + kthd1 * p1;
        kthd2 = jump_p * first_de2 + kthd2 * p1;
        kthd3 = jump_p * first_de3 + kthd3 * p1;
        
        arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
        arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
        arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
        
        
        jump_p=jump_p * p1;
        
        // Convert matrix to vector by row
        jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
        kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
        
        kthdvector += kthd1_vector(msmform(14)-1)*timeppvector;
        kthdvector1 += kthd2_vector(msmform(14)-1)*timeppvector;
        kthdvector2 += kthd3_vector(msmform(14)-1)*timeppvector;
      }
      timeMap[msmform(16)]= {kthdvector/kthpvector,kthdvector1/kthpvector, kthdvector2/kthpvector,
                             (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector,
                             (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector,
                             (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector,
                             (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector3/kthpvector,
                             (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector4/kthpvector,
                             (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector5/kthpvector
      };
      
      
      kthdvector_all += kthdvector/kthpvector;
      kthdvector1_all += kthdvector1/kthpvector;
      kthdvector2_all += kthdvector2/kthpvector;
      
      kthd2vector_all += (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector;
      kthd2vector1_all += (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector;
      kthd2vector2_all += (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector;
      kthd2vector3_all += (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector3/kthpvector;
      kthd2vector4_all += (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector4/kthpvector;
      kthd2vector5_all += (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector5/kthpvector;
      
    } else {
      
      
      kthdvector_all += timeMap[msmform(16)][0];
      kthdvector1_all += timeMap[msmform(16)][1];
      kthdvector2_all += timeMap[msmform(16)][2];
      
      kthd2vector_all += timeMap[msmform(16)][3];
      kthd2vector1_all += timeMap[msmform(16)][4];
      kthd2vector2_all += timeMap[msmform(16)][5];
      kthd2vector3_all += timeMap[msmform(16)][6];
      kthd2vector4_all += timeMap[msmform(16)][7];
      kthd2vector5_all += timeMap[msmform(16)][8];
      
    }
    
    
    
    //to_return_final[i]=kthpvector;
  }
  
  //arma::mat inform = {{kthd2vector_all, kthd2vector1_all, kthd2vector2_all},
  //                      {kthd2vector1_all, kthd2vector3_all, kthd2vector4_all},
  //                      {kthd2vector2_all, kthd2vector4_all, kthd2vector5_all}};
  
  //arma::mat score= {kthdvector_all,kthdvector1_all,kthdvector2_all};
  
  to_return_final ={kthdvector_all,kthdvector1_all,kthdvector2_all,kthd2vector_all,kthd2vector1_all,kthd2vector2_all,kthd2vector3_all,kthd2vector4_all,kthd2vector5_all};
  
  return to_return_final;
  
}

// [[Rcpp::export]]
std::vector<double> likelihood_part_no_map(arma::mat mat1) {
  
  
  int nrow_Mat = mat1.n_rows;
  std::vector<double>  to_return_final(9);
  
  double kthdvector_all  =0;
  double kthdvector1_all =0;
  double kthdvector2_all =0;
  
  double kthd2vector_all =0;
  double kthd2vector1_all=0;
  double kthd2vector2_all=0;
  double kthd2vector3_all=0;
  double kthd2vector4_all=0;
  double kthd2vector5_all=0;
  
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    arma::vec msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat jump_p=p1;
    
    arma::mat first_de1 = {{-1/msmform(12), 1/msmform(12), 0},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de2 = {{-1/msmform(12), 0, 1/msmform(12)},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de3 = {{0, 0, 0},
    {0, -1/msmform(12), 1/msmform(12)},
    {0, 0, 0}};
    
    
    
    arma::mat kthd1=first_de1;
    arma::mat kthd2=first_de2;
    arma::mat kthd3=first_de3;
    
    arma::mat kthd22 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_1 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_2 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_3 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_4 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_5 = arma::zeros<arma::mat>(3, 3);
    
    
    arma::mat jump_initial = arma::eye<arma::mat>(3, 3);
    
    
    // Create an identity matrix
    arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
    arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
    
    arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
    arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
    arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
    
    
    double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
      std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);
    
    
    
    double kthdvector = kthd1_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector1 = kthd2_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector2 = kthd3_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    
    double kthd2vector  = 0;
    double kthd2vector1 = 0;
    double kthd2vector2 = 0;
    double kthd2vector3 = 0;
    double kthd2vector4 = 0;
    double kthd2vector5 = 0;
    
    
    
    for (int j=0; j<20; j++) {
      //double timeppvector =std::exp(-msmform(12) * msmform(13)) * (std::pow(msmform(12) * msmform(13), j + 2)) / msm::factorial(j + 2);
      // numerical flow;
      double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
      
      kthd22 = kthd1 * first_de1 + kthd22 * p1 + kthd1 * first_de1 ;
      kthd22_1 = kthd2 * first_de1 + kthd22_1 * p1 + kthd1 * first_de2;
      kthd22_2 = kthd3 * first_de1 + kthd22_2 * p1 + kthd1 * first_de3;
      
      kthd22_3 = kthd2 * first_de2 + kthd22_3 * p1 + kthd2 * first_de2;
      kthd22_4 = kthd3 * first_de2 + kthd22_4 * p1 + kthd2 * first_de3;
      
      kthd22_5 = kthd3 * first_de3 + kthd22_5 * p1 + kthd3 * first_de3;
      
      arma::vec kthd22_vector = arma::vectorise(kthd22, 1).t();
      arma::vec kthd22_1_vector = arma::vectorise(kthd22_1,1).t();
      arma::vec kthd22_2_vector = arma::vectorise(kthd22_2, 1).t();
      arma::vec kthd22_3_vector = arma::vectorise(kthd22_3, 1).t();
      arma::vec kthd22_4_vector = arma::vectorise(kthd22_4, 1).t();
      arma::vec kthd22_5_vector = arma::vectorise(kthd22_5, 1).t();
      
      kthd2vector += kthd22_vector(msmform(14)-1) * timeppvector;
      kthd2vector1 += kthd22_1_vector(msmform(14)-1) * timeppvector;
      kthd2vector2 += kthd22_2_vector(msmform(14)-1) * timeppvector;
      kthd2vector3 += kthd22_3_vector(msmform(14)-1) * timeppvector;
      kthd2vector4 += kthd22_4_vector(msmform(14)-1) * timeppvector;
      kthd2vector5 += kthd22_5_vector(msmform(14)-1) * timeppvector;
      
      kthd1 = jump_p * first_de1 + kthd1 * p1;
      kthd2 = jump_p * first_de2 + kthd2 * p1;
      kthd3 = jump_p * first_de3 + kthd3 * p1;
      
      arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
      arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
      arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
      
      
      jump_p=jump_p * p1;
      
      // Convert matrix to vector by row
      jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
      kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
      
      kthdvector += kthd1_vector(msmform(14)-1)*timeppvector;
      kthdvector1 += kthd2_vector(msmform(14)-1)*timeppvector;
      kthdvector2 += kthd3_vector(msmform(14)-1)*timeppvector;
    }
    
    
    
    kthdvector_all += kthdvector/kthpvector;
    kthdvector1_all += kthdvector1/kthpvector;
    kthdvector2_all += kthdvector2/kthpvector;
    
    kthd2vector_all += (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector;
    kthd2vector1_all += (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector;
    kthd2vector2_all += (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector;
    kthd2vector3_all += (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector3/kthpvector;
    kthd2vector4_all += (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector4/kthpvector;
    kthd2vector5_all += (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector5/kthpvector;
    
    //to_return_final[i]=kthpvector;
  }
  
  //arma::mat inform = {{kthd2vector_all, kthd2vector1_all, kthd2vector2_all},
  //                      {kthd2vector1_all, kthd2vector3_all, kthd2vector4_all},
  //                      {kthd2vector2_all, kthd2vector4_all, kthd2vector5_all}};
  
  //arma::mat score= {kthdvector_all,kthdvector1_all,kthdvector2_all};
  
  to_return_final ={kthdvector_all,kthdvector1_all,kthdvector2_all,kthd2vector_all,kthd2vector1_all,kthd2vector2_all,kthd2vector3_all,kthd2vector4_all,kthd2vector5_all};
  
  return to_return_final;
  
}

// [[Rcpp::export]]
std::vector<double> likelihood_part2(arma::mat mat1) {
  
  
  int nrow_Mat = mat1.n_rows;
  std::vector<double>  to_return_final(9);
  
  double kthdvector_all  =0;
  double kthdvector1_all =0;
  double kthdvector2_all =0;
  double kthdvector3_all  =0;
  double kthdvector4_all =0;
  double kthdvector5_all =0;
  
  double kthd2vector_all =0;
  double kthd2vector1_all=0;
  double kthd2vector2_all=0;
  double kthd2vector3_all=0;
  double kthd2vector4_all=0;
  double kthd2vector5_all=0;
  
  double kthd2vector6_all=0;
  double kthd2vector7_all=0;
  double kthd2vector8_all=0;
  double kthd2vector9_all=0;
  double kthd2vector10_all=0;
  
  double kthd2vector11_all=0;
  double kthd2vector12_all=0;
  double kthd2vector13_all=0;
  double kthd2vector14_all=0;
  
  double kthd2vector15_all=0;
  double kthd2vector16_all=0;
  double kthd2vector17_all=0;
  
  double kthd2vector18_all=0;
  double kthd2vector19_all=0;
  
  double kthd2vector20_all=0;
  
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    arma::vec msmform= mat1.row(i).t() ;
    
    // Convert vector to matrix by row
    int numRows = 3;  // Number of rows in the resulting matrix
    int numCols = 3 ;  // Number of columns in the resulting matrix
    arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
    
    arma::mat jump_p=p1;
    
    arma::mat first_de1 = {{-1/msmform(12), 1/msmform(12), 0},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de2 = {{-1/msmform(12), 0, 1/msmform(12)},
    {0, 0, 0},
    {0, 0, 0}};
    
    arma::mat first_de3 = {{0, 0, 0},
    {1/msmform(12), -1/msmform(12), 0},
    {0, 0, 0}};
    
    arma::mat first_de4 = {{0, 0, 0},
    {0, -1/msmform(12), 1/msmform(12)},
    {0, 0, 0}};
    
    arma::mat first_de5 = {{0, 0, 0},
    {0, 0, 0},
    {1/msmform(12), 0, -1/msmform(12)}};
    
    arma::mat first_de6 = {{0, 0, 0},
    {0, 0, 0},
    {0, 1/msmform(12), -1/msmform(12)}};
    
    arma::mat kthd1=first_de1;
    arma::mat kthd2=first_de2;
    arma::mat kthd3=first_de3;
    arma::mat kthd4=first_de4;
    arma::mat kthd5=first_de5;
    arma::mat kthd6=first_de6;
    
    arma::mat kthd22 = arma::zeros<arma::mat>(3, 3);
    arma::mat kthd22_1 = kthd22 ;
    arma::mat kthd22_2 = kthd22 ;
    arma::mat kthd22_3 = kthd22 ;
    arma::mat kthd22_4 = kthd22 ;
    arma::mat kthd22_5 = kthd22 ;
    
    arma::mat kthd22_6 = kthd22 ;
    arma::mat kthd22_7 = kthd22 ;
    arma::mat kthd22_8 = kthd22 ;
    arma::mat kthd22_9 = kthd22 ;
    arma::mat kthd22_10 = kthd22 ;
    
    arma::mat kthd22_11 = kthd22 ;
    arma::mat kthd22_12 = kthd22 ;
    arma::mat kthd22_13 = kthd22 ;
    arma::mat kthd22_14 = kthd22 ;
    
    arma::mat kthd22_15 = kthd22 ;
    arma::mat kthd22_16 = kthd22 ;
    arma::mat kthd22_17 = kthd22 ;
    
    
    arma::mat kthd22_18 = kthd22 ;
    arma::mat kthd22_19 = kthd22 ;
    
    arma::mat kthd22_20 = kthd22 ;
    
    arma::mat jump_initial = arma::eye<arma::mat>(3, 3);
    
    // Create an identity matrix
    arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
    arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
    
    arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
    
    arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
    arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
    arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
    
    arma::vec kthd4_vector = arma::vectorise(kthd4, 1).t();
    arma::vec kthd5_vector = arma::vectorise(kthd5, 1).t();
    arma::vec kthd6_vector = arma::vectorise(kthd6, 1).t();
    
    
    double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
      std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);
    
    
    
    double kthdvector = kthd1_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector1 = kthd2_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector2 = kthd3_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector3 = kthd4_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector4 = kthd5_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    double kthdvector5 = kthd6_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
    
    double kthd2vector  = 0.0;
    double kthd2vector1 = 0.0;
    double kthd2vector2 = 0.0;
    double kthd2vector3 = 0.0;
    double kthd2vector4 = 0.0;
    double kthd2vector5 = 0.0;
    
    double kthd2vector6 = 0.0;
    double kthd2vector7 = 0.0;
    double kthd2vector8 = 0.0;
    double kthd2vector9 = 0.0;
    double kthd2vector10 = 0.0;
    
    double kthd2vector11 = 0.0;
    double kthd2vector12 = 0.0;
    double kthd2vector13 = 0.0;
    double kthd2vector14 = 0.0;
    
    double kthd2vector15 = 0.0;
    double kthd2vector16 = 0.0;
    double kthd2vector17 = 0.0;
    
    double kthd2vector18 = 0.0;
    double kthd2vector19 = 0.0;
    
    double kthd2vector20 = 0.0;
    
    for (int j=0; j<20; j++) {
      //double timeppvector =std::exp(-msmform(12) * msmform(13)) * (std::pow(msmform(12) * msmform(13), j + 2)) / msm::factorial(j + 2);
      // numerical flow;
      double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
      
      kthd22   = kthd1 * first_de1 + kthd22   * p1 + kthd1 * first_de1 ;
      kthd22_1 = kthd2 * first_de1 + kthd22_1 * p1 + kthd1 * first_de2;
      kthd22_2 = kthd3 * first_de1 + kthd22_2 * p1 + kthd1 * first_de3;
      kthd22_3 = kthd4 * first_de1 + kthd22_3 * p1 + kthd1 * first_de4;
      kthd22_4 = kthd5 * first_de1 + kthd22_4 * p1 + kthd1 * first_de5;
      kthd22_5 = kthd6 * first_de1 + kthd22_5 * p1 + kthd1 * first_de6;
      
      kthd22_6 = kthd2 * first_de2 + kthd22_6 * p1 + kthd2 * first_de2;
      kthd22_7 = kthd3 * first_de2 + kthd22_7 * p1 + kthd2 * first_de3;
      kthd22_8 = kthd4 * first_de2 + kthd22_8 * p1 + kthd2 * first_de4;
      kthd22_9 = kthd5 * first_de2 + kthd22_9 * p1 + kthd2 * first_de5;
      kthd22_10 = kthd6 * first_de2 + kthd22_10 * p1 + kthd2 * first_de6;
      
      kthd22_11 = kthd3 * first_de3 + kthd22_11 * p1 + kthd3 * first_de3;
      kthd22_12 = kthd4 * first_de3 + kthd22_12 * p1 + kthd3 * first_de4;
      kthd22_13 = kthd5 * first_de3 + kthd22_13 * p1 + kthd3 * first_de5;
      kthd22_14 = kthd6 * first_de3 + kthd22_14 * p1 + kthd3 * first_de6;
      
      kthd22_15 = kthd4 * first_de4 + kthd22_15 * p1 + kthd4 * first_de4;
      kthd22_16 = kthd5 * first_de4 + kthd22_16 * p1 + kthd4 * first_de5;
      kthd22_17 = kthd6 * first_de4 + kthd22_17 * p1 + kthd4 * first_de6;
      
      kthd22_18 = kthd5 * first_de5 + kthd22_18 * p1 + kthd5 * first_de5;
      kthd22_19 = kthd6 * first_de5 + kthd22_19 * p1 + kthd5 * first_de6;
      
      kthd22_20 = kthd6 * first_de6 + kthd22_20 * p1 + kthd6 * first_de6;
      
      
      arma::vec kthd22_vector = arma::vectorise(kthd22, 1).t();
      arma::vec kthd22_1_vector = arma::vectorise(kthd22_1,1).t();
      arma::vec kthd22_2_vector = arma::vectorise(kthd22_2, 1).t();
      arma::vec kthd22_3_vector = arma::vectorise(kthd22_3, 1).t();
      arma::vec kthd22_4_vector = arma::vectorise(kthd22_4, 1).t();
      arma::vec kthd22_5_vector = arma::vectorise(kthd22_5, 1).t();
      
      arma::vec kthd22_6_vector = arma::vectorise(kthd22_6,1).t();
      arma::vec kthd22_7_vector = arma::vectorise(kthd22_7, 1).t();
      arma::vec kthd22_8_vector = arma::vectorise(kthd22_8, 1).t();
      arma::vec kthd22_9_vector = arma::vectorise(kthd22_9, 1).t();
      arma::vec kthd22_10_vector = arma::vectorise(kthd22_10, 1).t();
      
      arma::vec kthd22_11_vector = arma::vectorise(kthd22_11, 1).t();
      arma::vec kthd22_12_vector = arma::vectorise(kthd22_12, 1).t();
      arma::vec kthd22_13_vector = arma::vectorise(kthd22_13, 1).t();
      arma::vec kthd22_14_vector = arma::vectorise(kthd22_14, 1).t();
      
      arma::vec kthd22_15_vector = arma::vectorise(kthd22_15, 1).t();
      arma::vec kthd22_16_vector = arma::vectorise(kthd22_16, 1).t();
      arma::vec kthd22_17_vector = arma::vectorise(kthd22_17, 1).t();
      
      arma::vec kthd22_18_vector = arma::vectorise(kthd22_18, 1).t();
      arma::vec kthd22_19_vector = arma::vectorise(kthd22_19, 1).t();
      
      arma::vec kthd22_20_vector = arma::vectorise(kthd22_20, 1).t();
      
      kthd2vector += kthd22_vector(msmform(14)-1) * timeppvector;
      kthd2vector1 += kthd22_1_vector(msmform(14)-1) * timeppvector;
      kthd2vector2 += kthd22_2_vector(msmform(14)-1) * timeppvector;
      kthd2vector3 += kthd22_3_vector(msmform(14)-1) * timeppvector;
      kthd2vector4 += kthd22_4_vector(msmform(14)-1) * timeppvector;
      kthd2vector5 += kthd22_5_vector(msmform(14)-1) * timeppvector;
      
      kthd2vector6 += kthd22_6_vector(msmform(14)-1) * timeppvector;
      kthd2vector7 += kthd22_7_vector(msmform(14)-1) * timeppvector;
      kthd2vector8 += kthd22_8_vector(msmform(14)-1) * timeppvector;
      kthd2vector9 += kthd22_9_vector(msmform(14)-1) * timeppvector;
      kthd2vector10 += kthd22_10_vector(msmform(14)-1) * timeppvector;
      
      kthd2vector11 += kthd22_11_vector(msmform(14)-1) * timeppvector;
      kthd2vector12 += kthd22_12_vector(msmform(14)-1) * timeppvector;
      kthd2vector13 += kthd22_13_vector(msmform(14)-1) * timeppvector;
      kthd2vector14 += kthd22_14_vector(msmform(14)-1) * timeppvector;
      
      kthd2vector15 += kthd22_15_vector(msmform(14)-1) * timeppvector;
      kthd2vector16 += kthd22_16_vector(msmform(14)-1) * timeppvector;
      kthd2vector17 += kthd22_17_vector(msmform(14)-1) * timeppvector;
      
      kthd2vector18 += kthd22_18_vector(msmform(14)-1) * timeppvector;
      kthd2vector19 += kthd22_19_vector(msmform(14)-1) * timeppvector;
      
      kthd2vector20 += kthd22_20_vector(msmform(14)-1) * timeppvector;
      
      kthd1 = jump_p * first_de1 + kthd1 * p1;
      kthd2 = jump_p * first_de2 + kthd2 * p1;
      kthd3 = jump_p * first_de3 + kthd3 * p1;
      
      kthd4 = jump_p * first_de4 + kthd4 * p1;
      kthd5 = jump_p * first_de5 + kthd5 * p1;
      kthd6 = jump_p * first_de6 + kthd6 * p1;
      
      arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
      arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
      arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
      
      arma::vec kthd4_vector = arma::vectorise(kthd4, 1).t();
      arma::vec kthd5_vector = arma::vectorise(kthd5, 1).t();
      arma::vec kthd6_vector = arma::vectorise(kthd6, 1).t();
      
      
      jump_p=jump_p * p1;
      
      // Convert matrix to vector by row
      jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
      kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
      
      kthdvector += kthd1_vector(msmform(14)-1)*timeppvector;
      kthdvector1 += kthd2_vector(msmform(14)-1)*timeppvector;
      kthdvector2 += kthd3_vector(msmform(14)-1)*timeppvector;
      
      kthdvector3 += kthd4_vector(msmform(14)-1)*timeppvector;
      kthdvector4 += kthd5_vector(msmform(14)-1)*timeppvector;
      kthdvector5 += kthd6_vector(msmform(14)-1)*timeppvector;
      
      
    }
    
    
    kthdvector_all += kthdvector/kthpvector;
    kthdvector1_all += kthdvector1/kthpvector;
    kthdvector2_all += kthdvector2/kthpvector;
    kthdvector3_all += kthdvector3/kthpvector;
    kthdvector4_all += kthdvector4/kthpvector;
    kthdvector5_all += kthdvector5/kthpvector;
    
    kthd2vector_all  += (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector;
    kthd2vector1_all += (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector;
    kthd2vector2_all += (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector;
    kthd2vector3_all += (-kthdvector*kthdvector3/(kthpvector*kthpvector))+kthd2vector3/kthpvector;
    kthd2vector4_all += (-kthdvector*kthdvector4/(kthpvector*kthpvector))+kthd2vector4/kthpvector;
    kthd2vector5_all += (-kthdvector*kthdvector5/(kthpvector*kthpvector))+kthd2vector5/kthpvector;
    
    kthd2vector6_all += (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector6/kthpvector;
    kthd2vector7_all += (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector7/kthpvector;
    kthd2vector8_all += (-kthdvector1*kthdvector3/(kthpvector*kthpvector))+kthd2vector8/kthpvector;
    kthd2vector9_all += (-kthdvector1*kthdvector4/(kthpvector*kthpvector))+kthd2vector9/kthpvector;
    kthd2vector10_all += (-kthdvector1*kthdvector5/(kthpvector*kthpvector))+kthd2vector10/kthpvector;
    
    kthd2vector11_all += (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector11/kthpvector;
    kthd2vector12_all += (-kthdvector2*kthdvector3/(kthpvector*kthpvector))+kthd2vector12/kthpvector;
    kthd2vector13_all += (-kthdvector2*kthdvector4/(kthpvector*kthpvector))+kthd2vector13/kthpvector;
    kthd2vector14_all += (-kthdvector2*kthdvector5/(kthpvector*kthpvector))+kthd2vector14/kthpvector;
    
    kthd2vector15_all += (-kthdvector3*kthdvector3/(kthpvector*kthpvector))+kthd2vector15/kthpvector;
    kthd2vector16_all += (-kthdvector3*kthdvector4/(kthpvector*kthpvector))+kthd2vector16/kthpvector;
    kthd2vector17_all += (-kthdvector3*kthdvector5/(kthpvector*kthpvector))+kthd2vector17/kthpvector;  
    
    kthd2vector18_all += (-kthdvector4*kthdvector4/(kthpvector*kthpvector))+kthd2vector18/kthpvector;
    kthd2vector19_all += (-kthdvector4*kthdvector5/(kthpvector*kthpvector))+kthd2vector19/kthpvector;  
    
    kthd2vector20_all += (-kthdvector5*kthdvector5/(kthpvector*kthpvector))+kthd2vector20/kthpvector; 
    //to_return_final[i]=kthpvector;
  }
  
  //arma::mat inform = {{kthd2vector_all, kthd2vector1_all, kthd2vector2_all},
  //                      {kthd2vector1_all, kthd2vector3_all, kthd2vector4_all},
  //                      {kthd2vector2_all, kthd2vector4_all, kthd2vector5_all}};
  
  //arma::mat score= {kthdvector_all,kthdvector1_all,kthdvector2_all};
  
  to_return_final ={kthdvector_all,kthdvector1_all,kthdvector2_all,kthdvector3_all,kthdvector4_all,kthdvector5_all,
                    kthd2vector_all,kthd2vector1_all,kthd2vector2_all,kthd2vector3_all,kthd2vector4_all,kthd2vector5_all,
                    kthd2vector6_all,kthd2vector7_all,kthd2vector8_all,kthd2vector9_all,kthd2vector10_all,kthd2vector11_all,
                    kthd2vector12_all,kthd2vector13_all,kthd2vector14_all,kthd2vector15_all,kthd2vector16_all,kthd2vector17_all,
                    kthd2vector18_all,kthd2vector19_all,kthd2vector20_all};
  
  return to_return_final;
  
}

// [[Rcpp::export]]
std::vector<double> likelihood_part2_map(arma::mat mat1) {
  
  std::map<int, std::vector<double>> timeMap;
  int nrow_Mat = mat1.n_rows;
  std::vector<double>  to_return_final(9);
  
  double kthdvector_all  =0;
  double kthdvector1_all =0;
  double kthdvector2_all =0;
  double kthdvector3_all  =0;
  double kthdvector4_all =0;
  double kthdvector5_all =0;
  
  double kthd2vector_all =0;
  double kthd2vector1_all=0;
  double kthd2vector2_all=0;
  double kthd2vector3_all=0;
  double kthd2vector4_all=0;
  double kthd2vector5_all=0;
  
  double kthd2vector6_all=0;
  double kthd2vector7_all=0;
  double kthd2vector8_all=0;
  double kthd2vector9_all=0;
  double kthd2vector10_all=0;
  
  double kthd2vector11_all=0;
  double kthd2vector12_all=0;
  double kthd2vector13_all=0;
  double kthd2vector14_all=0;
  
  double kthd2vector15_all=0;
  double kthd2vector16_all=0;
  double kthd2vector17_all=0;
  
  double kthd2vector18_all=0;
  double kthd2vector19_all=0;
  
  double kthd2vector20_all=0;
  
  
  for ( int i=0; i<nrow_Mat; i++ ) {
    arma::vec msmform= mat1.row(i).t() ;
    
    if (timeMap.count(msmform(16)) == 0) { 
      
      // Convert vector to matrix by row
      int numRows = 3;  // Number of rows in the resulting matrix
      int numCols = 3 ;  // Number of columns in the resulting matrix
      arma::mat p1 = arma::reshape(msmform.subvec(0, 8), numRows, numCols).t();
      
      arma::mat jump_p=p1;
      
      arma::mat first_de1 = {{-1/msmform(12), 1/msmform(12), 0},
      {0, 0, 0},
      {0, 0, 0}};
      
      arma::mat first_de2 = {{-1/msmform(12), 0, 1/msmform(12)},
      {0, 0, 0},
      {0, 0, 0}};
      
      arma::mat first_de3 = {{0, 0, 0},
      {1/msmform(12), -1/msmform(12), 0},
      {0, 0, 0}};
      
      arma::mat first_de4 = {{0, 0, 0},
      {0, -1/msmform(12), 1/msmform(12)},
      {0, 0, 0}};
      
      arma::mat first_de5 = {{0, 0, 0},
      {0, 0, 0},
      {1/msmform(12), 0, -1/msmform(12)}};
      
      arma::mat first_de6 = {{0, 0, 0},
      {0, 0, 0},
      {0, 1/msmform(12), -1/msmform(12)}};
      
      arma::mat kthd1=first_de1;
      arma::mat kthd2=first_de2;
      arma::mat kthd3=first_de3;
      arma::mat kthd4=first_de4;
      arma::mat kthd5=first_de5;
      arma::mat kthd6=first_de6;
      
      arma::mat kthd22 = arma::zeros<arma::mat>(3, 3);
      arma::mat kthd22_1 = kthd22 ;
      arma::mat kthd22_2 = kthd22 ;
      arma::mat kthd22_3 = kthd22 ;
      arma::mat kthd22_4 = kthd22 ;
      arma::mat kthd22_5 = kthd22 ;
      
      arma::mat kthd22_6 = kthd22 ;
      arma::mat kthd22_7 = kthd22 ;
      arma::mat kthd22_8 = kthd22 ;
      arma::mat kthd22_9 = kthd22 ;
      arma::mat kthd22_10 = kthd22 ;
      
      arma::mat kthd22_11 = kthd22 ;
      arma::mat kthd22_12 = kthd22 ;
      arma::mat kthd22_13 = kthd22 ;
      arma::mat kthd22_14 = kthd22 ;
      
      arma::mat kthd22_15 = kthd22 ;
      arma::mat kthd22_16 = kthd22 ;
      arma::mat kthd22_17 = kthd22 ;
      
      
      arma::mat kthd22_18 = kthd22 ;
      arma::mat kthd22_19 = kthd22 ;
      
      arma::mat kthd22_20 = kthd22 ;
      
      arma::mat jump_initial = arma::eye<arma::mat>(3, 3);
      
      // Create an identity matrix
      arma::mat identity_mat1 = arma::eye<arma::mat>(3, 3);
      arma::vec identity_mat1_vector = arma::vectorise(identity_mat1, 1).t();
      
      arma::vec jump_p_vector = arma::vectorise(jump_p, 1).t();
      
      arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
      arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
      arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
      
      arma::vec kthd4_vector = arma::vectorise(kthd4, 1).t();
      arma::vec kthd5_vector = arma::vectorise(kthd5, 1).t();
      arma::vec kthd6_vector = arma::vectorise(kthd6, 1).t();
      
      
      double kthpvector = std::exp(-msmform(12) * msmform(13)) * identity_mat1_vector(msmform(14)-1) +
        std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13) * jump_p_vector(msmform(14)-1);
      
      
      
      double kthdvector = kthd1_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      double kthdvector1 = kthd2_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      double kthdvector2 = kthd3_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      double kthdvector3 = kthd4_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      double kthdvector4 = kthd5_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      double kthdvector5 = kthd6_vector(msmform(14)-1) * std::exp(-msmform(12) * msmform(13)) * msmform(12) * msmform(13);
      
      double kthd2vector  = 0.0;
      double kthd2vector1 = 0.0;
      double kthd2vector2 = 0.0;
      double kthd2vector3 = 0.0;
      double kthd2vector4 = 0.0;
      double kthd2vector5 = 0.0;
      
      double kthd2vector6 = 0.0;
      double kthd2vector7 = 0.0;
      double kthd2vector8 = 0.0;
      double kthd2vector9 = 0.0;
      double kthd2vector10 = 0.0;
      
      double kthd2vector11 = 0.0;
      double kthd2vector12 = 0.0;
      double kthd2vector13 = 0.0;
      double kthd2vector14 = 0.0;
      
      double kthd2vector15 = 0.0;
      double kthd2vector16 = 0.0;
      double kthd2vector17 = 0.0;
      
      double kthd2vector18 = 0.0;
      double kthd2vector19 = 0.0;
      
      double kthd2vector20 = 0.0;
      
      for (int j=0; j<20; j++) {
        //double timeppvector =std::exp(-msmform(12) * msmform(13)) * (std::pow(msmform(12) * msmform(13), j + 2)) / msm::factorial(j + 2);
        // numerical flow;
        double timeppvector =std::exp((-msmform(12) * msmform(13)) + (j + 2) * log(msmform(12) * msmform(13)) - (msm::logfactorial(j+2)));
        
        kthd22   = kthd1 * first_de1 + kthd22   * p1 + kthd1 * first_de1 ;
        kthd22_1 = kthd2 * first_de1 + kthd22_1 * p1 + kthd1 * first_de2;
        kthd22_2 = kthd3 * first_de1 + kthd22_2 * p1 + kthd1 * first_de3;
        kthd22_3 = kthd4 * first_de1 + kthd22_3 * p1 + kthd1 * first_de4;
        kthd22_4 = kthd5 * first_de1 + kthd22_4 * p1 + kthd1 * first_de5;
        kthd22_5 = kthd6 * first_de1 + kthd22_5 * p1 + kthd1 * first_de6;
        
        kthd22_6 = kthd2 * first_de2 + kthd22_6 * p1 + kthd2 * first_de2;
        kthd22_7 = kthd3 * first_de2 + kthd22_7 * p1 + kthd2 * first_de3;
        kthd22_8 = kthd4 * first_de2 + kthd22_8 * p1 + kthd2 * first_de4;
        kthd22_9 = kthd5 * first_de2 + kthd22_9 * p1 + kthd2 * first_de5;
        kthd22_10 = kthd6 * first_de2 + kthd22_10 * p1 + kthd2 * first_de6;
        
        kthd22_11 = kthd3 * first_de3 + kthd22_11 * p1 + kthd3 * first_de3;
        kthd22_12 = kthd4 * first_de3 + kthd22_12 * p1 + kthd3 * first_de4;
        kthd22_13 = kthd5 * first_de3 + kthd22_13 * p1 + kthd3 * first_de5;
        kthd22_14 = kthd6 * first_de3 + kthd22_14 * p1 + kthd3 * first_de6;
        
        kthd22_15 = kthd4 * first_de4 + kthd22_15 * p1 + kthd4 * first_de4;
        kthd22_16 = kthd5 * first_de4 + kthd22_16 * p1 + kthd4 * first_de5;
        kthd22_17 = kthd6 * first_de4 + kthd22_17 * p1 + kthd4 * first_de6;
        
        kthd22_18 = kthd5 * first_de5 + kthd22_18 * p1 + kthd5 * first_de5;
        kthd22_19 = kthd6 * first_de5 + kthd22_19 * p1 + kthd5 * first_de6;
        
        kthd22_20 = kthd6 * first_de6 + kthd22_20 * p1 + kthd6 * first_de6;
        
        
        arma::vec kthd22_vector = arma::vectorise(kthd22, 1).t();
        arma::vec kthd22_1_vector = arma::vectorise(kthd22_1,1).t();
        arma::vec kthd22_2_vector = arma::vectorise(kthd22_2, 1).t();
        arma::vec kthd22_3_vector = arma::vectorise(kthd22_3, 1).t();
        arma::vec kthd22_4_vector = arma::vectorise(kthd22_4, 1).t();
        arma::vec kthd22_5_vector = arma::vectorise(kthd22_5, 1).t();
        
        arma::vec kthd22_6_vector = arma::vectorise(kthd22_6,1).t();
        arma::vec kthd22_7_vector = arma::vectorise(kthd22_7, 1).t();
        arma::vec kthd22_8_vector = arma::vectorise(kthd22_8, 1).t();
        arma::vec kthd22_9_vector = arma::vectorise(kthd22_9, 1).t();
        arma::vec kthd22_10_vector = arma::vectorise(kthd22_10, 1).t();
        
        arma::vec kthd22_11_vector = arma::vectorise(kthd22_11, 1).t();
        arma::vec kthd22_12_vector = arma::vectorise(kthd22_12, 1).t();
        arma::vec kthd22_13_vector = arma::vectorise(kthd22_13, 1).t();
        arma::vec kthd22_14_vector = arma::vectorise(kthd22_14, 1).t();
        
        arma::vec kthd22_15_vector = arma::vectorise(kthd22_15, 1).t();
        arma::vec kthd22_16_vector = arma::vectorise(kthd22_16, 1).t();
        arma::vec kthd22_17_vector = arma::vectorise(kthd22_17, 1).t();
        
        arma::vec kthd22_18_vector = arma::vectorise(kthd22_18, 1).t();
        arma::vec kthd22_19_vector = arma::vectorise(kthd22_19, 1).t();
        
        arma::vec kthd22_20_vector = arma::vectorise(kthd22_20, 1).t();
        
        kthd2vector += kthd22_vector(msmform(14)-1) * timeppvector;
        kthd2vector1 += kthd22_1_vector(msmform(14)-1) * timeppvector;
        kthd2vector2 += kthd22_2_vector(msmform(14)-1) * timeppvector;
        kthd2vector3 += kthd22_3_vector(msmform(14)-1) * timeppvector;
        kthd2vector4 += kthd22_4_vector(msmform(14)-1) * timeppvector;
        kthd2vector5 += kthd22_5_vector(msmform(14)-1) * timeppvector;
        
        kthd2vector6 += kthd22_6_vector(msmform(14)-1) * timeppvector;
        kthd2vector7 += kthd22_7_vector(msmform(14)-1) * timeppvector;
        kthd2vector8 += kthd22_8_vector(msmform(14)-1) * timeppvector;
        kthd2vector9 += kthd22_9_vector(msmform(14)-1) * timeppvector;
        kthd2vector10 += kthd22_10_vector(msmform(14)-1) * timeppvector;
        
        kthd2vector11 += kthd22_11_vector(msmform(14)-1) * timeppvector;
        kthd2vector12 += kthd22_12_vector(msmform(14)-1) * timeppvector;
        kthd2vector13 += kthd22_13_vector(msmform(14)-1) * timeppvector;
        kthd2vector14 += kthd22_14_vector(msmform(14)-1) * timeppvector;
        
        kthd2vector15 += kthd22_15_vector(msmform(14)-1) * timeppvector;
        kthd2vector16 += kthd22_16_vector(msmform(14)-1) * timeppvector;
        kthd2vector17 += kthd22_17_vector(msmform(14)-1) * timeppvector;
        
        kthd2vector18 += kthd22_18_vector(msmform(14)-1) * timeppvector;
        kthd2vector19 += kthd22_19_vector(msmform(14)-1) * timeppvector;
        
        kthd2vector20 += kthd22_20_vector(msmform(14)-1) * timeppvector;
        
        kthd1 = jump_p * first_de1 + kthd1 * p1;
        kthd2 = jump_p * first_de2 + kthd2 * p1;
        kthd3 = jump_p * first_de3 + kthd3 * p1;
        
        kthd4 = jump_p * first_de4 + kthd4 * p1;
        kthd5 = jump_p * first_de5 + kthd5 * p1;
        kthd6 = jump_p * first_de6 + kthd6 * p1;
        
        arma::vec kthd1_vector = arma::vectorise(kthd1, 1).t();
        arma::vec kthd2_vector = arma::vectorise(kthd2, 1).t();
        arma::vec kthd3_vector = arma::vectorise(kthd3, 1).t();
        
        arma::vec kthd4_vector = arma::vectorise(kthd4, 1).t();
        arma::vec kthd5_vector = arma::vectorise(kthd5, 1).t();
        arma::vec kthd6_vector = arma::vectorise(kthd6, 1).t();
        
        
        jump_p=jump_p * p1;
        
        // Convert matrix to vector by row
        jump_p_vector = arma::vectorise(jump_p, 1).t(); // 1 indicates row-wise vectorisation
        kthpvector += jump_p_vector(msmform(14)-1)*timeppvector;
        
        kthdvector += kthd1_vector(msmform(14)-1)*timeppvector;
        kthdvector1 += kthd2_vector(msmform(14)-1)*timeppvector;
        kthdvector2 += kthd3_vector(msmform(14)-1)*timeppvector;
        
        kthdvector3 += kthd4_vector(msmform(14)-1)*timeppvector;
        kthdvector4 += kthd5_vector(msmform(14)-1)*timeppvector;
        kthdvector5 += kthd6_vector(msmform(14)-1)*timeppvector;
        
      } 
      
      
      kthdvector_all += kthdvector/kthpvector;
      kthdvector1_all += kthdvector1/kthpvector;
      kthdvector2_all += kthdvector2/kthpvector;
      kthdvector3_all += kthdvector3/kthpvector;
      kthdvector4_all += kthdvector4/kthpvector;
      kthdvector5_all += kthdvector5/kthpvector;
      
      kthd2vector_all  += (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector;
      kthd2vector1_all += (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector;
      kthd2vector2_all += (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector;
      kthd2vector3_all += (-kthdvector*kthdvector3/(kthpvector*kthpvector))+kthd2vector3/kthpvector;
      kthd2vector4_all += (-kthdvector*kthdvector4/(kthpvector*kthpvector))+kthd2vector4/kthpvector;
      kthd2vector5_all += (-kthdvector*kthdvector5/(kthpvector*kthpvector))+kthd2vector5/kthpvector;
      
      kthd2vector6_all += (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector6/kthpvector;
      kthd2vector7_all += (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector7/kthpvector;
      kthd2vector8_all += (-kthdvector1*kthdvector3/(kthpvector*kthpvector))+kthd2vector8/kthpvector;
      kthd2vector9_all += (-kthdvector1*kthdvector4/(kthpvector*kthpvector))+kthd2vector9/kthpvector;
      kthd2vector10_all += (-kthdvector1*kthdvector5/(kthpvector*kthpvector))+kthd2vector10/kthpvector;
      
      kthd2vector11_all += (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector11/kthpvector;
      kthd2vector12_all += (-kthdvector2*kthdvector3/(kthpvector*kthpvector))+kthd2vector12/kthpvector;
      kthd2vector13_all += (-kthdvector2*kthdvector4/(kthpvector*kthpvector))+kthd2vector13/kthpvector;
      kthd2vector14_all += (-kthdvector2*kthdvector5/(kthpvector*kthpvector))+kthd2vector14/kthpvector;
      
      kthd2vector15_all += (-kthdvector3*kthdvector3/(kthpvector*kthpvector))+kthd2vector15/kthpvector;
      kthd2vector16_all += (-kthdvector3*kthdvector4/(kthpvector*kthpvector))+kthd2vector16/kthpvector;
      kthd2vector17_all += (-kthdvector3*kthdvector5/(kthpvector*kthpvector))+kthd2vector17/kthpvector;  
      
      kthd2vector18_all += (-kthdvector4*kthdvector4/(kthpvector*kthpvector))+kthd2vector18/kthpvector;
      kthd2vector19_all += (-kthdvector4*kthdvector5/(kthpvector*kthpvector))+kthd2vector19/kthpvector;  
      
      kthd2vector20_all += (-kthdvector5*kthdvector5/(kthpvector*kthpvector))+kthd2vector20/kthpvector; 
      
      timeMap[msmform(16)]={
        
        kthdvector/kthpvector,
        kthdvector1/kthpvector,
        kthdvector2/kthpvector,
        kthdvector3/kthpvector,
        kthdvector4/kthpvector,
        kthdvector5/kthpvector,
        
        (-kthdvector*kthdvector/(kthpvector*kthpvector))+kthd2vector/kthpvector,
        (-kthdvector*kthdvector1/(kthpvector*kthpvector))+kthd2vector1/kthpvector,
        (-kthdvector*kthdvector2/(kthpvector*kthpvector))+kthd2vector2/kthpvector,
        (-kthdvector*kthdvector3/(kthpvector*kthpvector))+kthd2vector3/kthpvector,
        (-kthdvector*kthdvector4/(kthpvector*kthpvector))+kthd2vector4/kthpvector,
        (-kthdvector*kthdvector5/(kthpvector*kthpvector))+kthd2vector5/kthpvector,
        
        (-kthdvector1*kthdvector1/(kthpvector*kthpvector))+kthd2vector6/kthpvector,
        (-kthdvector1*kthdvector2/(kthpvector*kthpvector))+kthd2vector7/kthpvector,
        (-kthdvector1*kthdvector3/(kthpvector*kthpvector))+kthd2vector8/kthpvector,
        (-kthdvector1*kthdvector4/(kthpvector*kthpvector))+kthd2vector9/kthpvector,
        (-kthdvector1*kthdvector5/(kthpvector*kthpvector))+kthd2vector10/kthpvector,
        
        (-kthdvector2*kthdvector2/(kthpvector*kthpvector))+kthd2vector11/kthpvector,
        (-kthdvector2*kthdvector3/(kthpvector*kthpvector))+kthd2vector12/kthpvector,
        (-kthdvector2*kthdvector4/(kthpvector*kthpvector))+kthd2vector13/kthpvector,
        (-kthdvector2*kthdvector5/(kthpvector*kthpvector))+kthd2vector14/kthpvector,
        
        (-kthdvector3*kthdvector3/(kthpvector*kthpvector))+kthd2vector15/kthpvector,
        (-kthdvector3*kthdvector4/(kthpvector*kthpvector))+kthd2vector16/kthpvector,
        (-kthdvector3*kthdvector5/(kthpvector*kthpvector))+kthd2vector17/kthpvector,  
        
        (-kthdvector4*kthdvector4/(kthpvector*kthpvector))+kthd2vector18/kthpvector,
        (-kthdvector4*kthdvector5/(kthpvector*kthpvector))+kthd2vector19/kthpvector,  
        
        (-kthdvector5*kthdvector5/(kthpvector*kthpvector))+kthd2vector20/kthpvector, 
        
      };
      
    }else {
      
      
      kthdvector_all += timeMap[msmform(16)][0];
      kthdvector1_all += timeMap[msmform(16)][1];
      kthdvector2_all += timeMap[msmform(16)][2];
      kthdvector3_all += timeMap[msmform(16)][3];
      kthdvector4_all += timeMap[msmform(16)][4];
      kthdvector5_all += timeMap[msmform(16)][5];
      
      kthd2vector_all  += timeMap[msmform(16)][6];
      kthd2vector1_all += timeMap[msmform(16)][7];
      kthd2vector2_all += timeMap[msmform(16)][8];
      kthd2vector3_all += timeMap[msmform(16)][9];
      kthd2vector4_all += timeMap[msmform(16)][10];
      kthd2vector5_all += timeMap[msmform(16)][11];
      
      kthd2vector6_all += timeMap[msmform(16)][12];
      kthd2vector7_all += timeMap[msmform(16)][13];
      kthd2vector8_all += timeMap[msmform(16)][14];
      kthd2vector9_all += timeMap[msmform(16)][15];
      kthd2vector10_all += timeMap[msmform(16)][16];
      
      kthd2vector11_all += timeMap[msmform(16)][17];
      kthd2vector12_all += timeMap[msmform(16)][18];
      kthd2vector13_all += timeMap[msmform(16)][19];
      kthd2vector14_all += timeMap[msmform(16)][20];
      
      kthd2vector15_all += timeMap[msmform(16)][21];
      kthd2vector16_all += timeMap[msmform(16)][22];
      kthd2vector17_all +=timeMap[msmform(16)][23];  
      
      kthd2vector18_all += timeMap[msmform(16)][24];
      kthd2vector19_all += timeMap[msmform(16)][25];  
      
      kthd2vector20_all += timeMap[msmform(16)][26]; 
      
    }
    
  }
  
  to_return_final ={kthdvector_all,kthdvector1_all,kthdvector2_all,kthdvector3_all,kthdvector4_all,kthdvector5_all,
                    kthd2vector_all,kthd2vector1_all,kthd2vector2_all,kthd2vector3_all,kthd2vector4_all,kthd2vector5_all,
                    kthd2vector6_all,kthd2vector7_all,kthd2vector8_all,kthd2vector9_all,kthd2vector10_all,kthd2vector11_all,
                    kthd2vector12_all,kthd2vector13_all,kthd2vector14_all,kthd2vector15_all,kthd2vector16_all,kthd2vector17_all,
                    kthd2vector18_all,kthd2vector19_all,kthd2vector20_all};
  
  return to_return_final;
  
}






