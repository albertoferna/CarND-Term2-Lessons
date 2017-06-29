#include <iostream>
#include "ukf.h"

UKF::UKF() {
  //TODO Auto-generated constructor stub
  Init();
}

UKF::~UKF() {
  //TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
     Xsig_aug <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  double delta_t = 0.1; //time diff in sec
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  Xsig_pred.fill(0.0);
  VectorXd aug_state_k(7);
  VectorXd noise_k(5);
  for (int i=0; i< 2 * n_aug + 1; i++) {
      aug_state_k = Xsig_aug.col(i);
      noise_k.fill(0.0);
      noise_k(0) = 0.5 * delta_t * delta_t * cos(aug_state_k(3)) * aug_state_k(5);
      noise_k(1) = 0.5 * delta_t * delta_t * sin(aug_state_k(3)) * aug_state_k(5);
      noise_k(2) = delta_t * aug_state_k(5);
      noise_k(3) = 0.5 * delta_t * delta_t * aug_state_k(6);
      noise_k(4) = delta_t * aug_state_k(6);
      if (aug_state_k(4) < 0.00001) {
          Xsig_pred.col(i) = aug_state_k.head(5) + noise_k;
          Xsig_pred(0, i) += aug_state_k(2) * cos(aug_state_k(3)) * delta_t;
          Xsig_pred(1, i) += aug_state_k(2) * sin(aug_state_k(3)) * delta_t;
          Xsig_pred(3, i) += aug_state_k(4) * delta_t;
      } else {
          Xsig_pred.col(i) = aug_state_k.head(5) + noise_k;
          Xsig_pred(0, i) += aug_state_k(2)/aug_state_k(4) * (sin(aug_state_k(3) + aug_state_k(4) * delta_t) - sin(aug_state_k(3)));
          Xsig_pred(1, i) += aug_state_k(2)/aug_state_k(4) * (-cos(aug_state_k(3) + aug_state_k(4) * delta_t) + cos(aug_state_k(3)));
          Xsig_pred(3, i) += aug_state_k(4) * delta_t;
      }
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}
