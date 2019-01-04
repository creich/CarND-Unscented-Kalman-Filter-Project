#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 3;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  P_.setIdentity();
  //std::cout << "P" << std::endl;
  //std::cout << P_ << std::endl;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

   if( !use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER ) return;
   if( !use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR ) return;

   // init
   if( false == is_initialized_ ) {
      // TODO

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
         // TODO: Convert radar from polar to cartesian coordinates
         //         and initialize state.
         float roh     = meas_package.raw_measurements_(0);
         float phi     = meas_package.raw_measurements_(1);
         float roh_dot = meas_package.raw_measurements_(2);

         x_ << roh * cos(phi),
               roh * sin(phi),
               roh_dot * cos(phi),
               roh_dot * sin(phi),
               0.0;
       }
       else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
         // TODO: Initialize state.
         // set the state with the initial location and zero velocity
         x_ << meas_package.raw_measurements_(0),
               meas_package.raw_measurements_(1),
               0.0,
               0.0,
               0.0;
      }

      time_us_ = meas_package.timestamp_;

      // done initializing, no need to predict or update
      is_initialized_ = true;

      return;
   }

   double delta_t = (time_us_ - meas_package.timestamp_) / 1000000;

   // predict
//   std::cout << "###################" << std::endl;
//   std::cout << "x_ before" << std::endl;
//   std::cout << x_ << std::endl;

   Prediction(delta_t);

//   std::cout << "x_ after" << std::endl;
//   std::cout << x_ << std::endl;

   // update
   switch (meas_package.sensor_type_) {
      case MeasurementPackage::LASER:
         UpdateLidar(meas_package);
         break;
      case MeasurementPackage::RADAR:
         UpdateRadar(meas_package);
         break;
   }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //##############  generate sigma points  ###############

  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  // calculate sigma points
  // first column
  Xsig.col(0) = x_;

  MatrixXd A = (lambda_ + n_x_) * P_;
  A = A.llt().matrixL();

  // set sigma points as columns of matrix Xsig
  for(int i=0; i<n_x_; i++) {
      // second to n_x+1
      Xsig.col(1+i) = x_ + A.col(i);
      // n_x+2 to 2*n_x+1
      Xsig.col(n_x_+1+i) = x_ - A.col(i);
  }

  //##############  augment sigma points  ###############

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  int n_diff = n_aug_ - n_x_;
 
  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug.tail(n_diff) << 0, 0;

  // create augmented covariance matrix
     // process noise covariance matrix
  MatrixXd Q = MatrixXd(n_diff, n_diff);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
       
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_diff, n_diff) = Q;

  // first column
  Xsig_aug.col(0) = x_aug;
  
  MatrixXd B = (lambda_ + n_aug_) * P_aug;
  // create square root matrix
  B = B.llt().matrixL();

  // create augmented sigma points
  for(int i=0; i<n_aug_; i++) {
      // second to n_aug+1
      Xsig_aug.col(1+i) = x_aug + B.col(i);
      // n_x+2 to 2*n_aug+1
      Xsig_aug.col(n_aug_+1+i) = x_aug - B.col(i);
  }

  //##############  predict sigma points  ###############

  // create matrix for predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  for(int i=0; i < 2 * n_aug_ + 1; i++) {
      VectorXd x_k(5);  // x_k
      VectorXd x_k1(5); // summand
      VectorXd v_k(5);  // process noise

      // get values for calculation --> x_k
      x_k = Xsig_aug.col(i).topRows(5);
      // for the ease of calculation...
      float px = x_k(0);
      float py = x_k(1);
      float v  = x_k(2);
      float psi = x_k(3);
      float psi_dot = x_k(4);
      // noise values
      float v_a = Xsig_aug.col(i)(5);
      float v_psi_dotdot = Xsig_aug.col(i)(6);

      // avoid division by zero
      // TODO check if this is sufficient (due to float values not exactly 0)
      if(fabs(psi_dot) < 0.0001) {
          // summand (part 1)
          x_k1(0) = v * cos(psi) * delta_t;
          x_k1(1) = v * sin(psi) * delta_t;
      } else {
          // summand (part 1)
          x_k1(0) = v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi));
          x_k1(1) = v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi));
      }
      // summand (part 2)
      x_k1(2) = 0;
      x_k1(3) = psi_dot * delta_t;
      x_k1(4) = 0;

      // noise vector
      v_k(0) = 0.5 * delta_t*delta_t * cos(psi) * v_a;
      v_k(1) = 0.5 * delta_t*delta_t * sin(psi) * v_a;
      v_k(2) = delta_t * v_a;
      v_k(3) = 0.5 * delta_t*delta_t * v_psi_dotdot;
      v_k(4) = delta_t * v_psi_dotdot;

      // write predicted sigma points into right column
      Xsig_pred_.col(i) = x_k + x_k1 + v_k;
  }

  //##############  predict mean and covariance  ###############

  // create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  // set weights
  int n_a = 2 * n_aug_ + 1;
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i=1; i<n_a; i++) {
      weights_(i) = 1. / (2. * (lambda_ + n_aug_));
  }

  // predict state mean
  x_.fill(0.0);
  for(int i=0; i<n_a; i++) {
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  P_.fill(0.0);
  for(int i=0; i<n_a; i++) {
      VectorXd current_vec(n_x_);
      current_vec = Xsig_pred_.col(i) - x_;
      // TODO check angle normalization from tutors sollution
      P_ = P_ + ( weights_(i) * (current_vec * current_vec.transpose()) );
  }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_x_);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_x_,n_x_);

  // calculate mean predicted measurement
  z_pred.fill(0.0);

  for (int i=0; i < 2*n_aug_+1; ++i) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S

  S.fill(0.0);
  for(int i=0; i<2 * n_aug_ + 1; i++) {
      VectorXd current_vec(n_x_);
      current_vec = Zsig.col(i) - z_pred;

      // angle normalization
      while (current_vec(1)> M_PI) current_vec(1)-=2.*M_PI;
      while (current_vec(1)<-M_PI) current_vec(1)+=2.*M_PI;

      S = S + ( weights_(i) * current_vec * current_vec.transpose() );
  }

  MatrixXd R = MatrixXd(n_x_,n_x_);
  R.fill(0.0);
  R(0, 0) = std_laspx_*std_laspx_;
  R(1, 1) = std_laspy_*std_laspy_;
  R(2, 2) = std_a_*std_a_;
  R(3, 3) = std_yawdd_*std_yawdd_;
//  R << std_laspx_*std_laspx_, 0, 0, 0, 0
//       0, std_laspy_*std_laspy_, 0, 0, 0,
//       0, 0, std_a_*std_a_, 0, 0,
//       0, 0, 0, std_yawdd_*std_yawdd_, 0,
//       0, 0, 0, 0, 0;

std::cout << "S" << std::endl;
std::cout << S << std::endl;
std::cout << "R" << std::endl;
std::cout << R << std::endl;

  S = S + R;

  //############# use incoming data ####################

  // create example vector for incoming lidar measurement
  VectorXd z = VectorXd(n_x_);
  z.fill(0.0);
  z.topRows(2) = meas_package.raw_measurements_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_x_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i<2 * n_aug_ + 1; i++) {

std::cout << "Zsig" << std::endl;
std::cout << Zsig << std::endl;
std::cout << "z_pred" << std::endl;
std::cout << z_pred << std::endl;

      // residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

std::cout << "Xsig_pred_" << std::endl;
std::cout << Xsig_pred_ << std::endl;
std::cout << "z" << std::endl;
std::cout << z << std::endl;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - z;
      // angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
/*
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
*/
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  //################## predict radar data #######################

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  for(int i=0; i<2 * n_aug_ + 1; i++) {
      // pick data from Xsig_pred
      double px      = Xsig_pred_.col(i)(0);
      double py      = Xsig_pred_.col(i)(1);
      double v       = Xsig_pred_.col(i)(2);
      double psi     = Xsig_pred_.col(i)(3);
      double psi_dot = Xsig_pred_.col(i)(4);

      // transform dadta and store it into Zsig
      Zsig.col(i) << sqrt(px*px + py*py),
                     atan2(py, px),
                     (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px*px + py*py);
  }

  // calculate mean predicted measurement
  // following stuff (for some unknown reason) gives a slightly different result
  // z_pred = Zsig.rowwise().mean();
  z_pred.fill(0.0);
/*
  std::cout << "z_pred" << std::endl;
  std::cout << z_pred << std::endl;
  std::cout << "weights_" << std::endl;
  std::cout << weights_ << std::endl;
  std::cout << "Zsig" << std::endl;
  std::cout << Zsig << std::endl;
*/
  for (int i=0; i < 2*n_aug_+1; ++i) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S

  S.fill(0.0);
  for(int i=0; i<2 * n_aug_ + 1; i++) {
      VectorXd current_vec(n_x_);
      current_vec = Zsig.col(i) - z_pred;

      // angle normalization
      while (current_vec(1)> M_PI) current_vec(1)-=2.*M_PI;
      while (current_vec(1)<-M_PI) current_vec(1)+=2.*M_PI;

      S = S + ( weights_(i) * current_vec * current_vec.transpose() );
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S = S + R;

  //############# use incoming data ####################

  // create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i<2 * n_aug_ + 1; i++) {
      // TODO introduce angle normalization from tutors solution
      //Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();

      // residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      // angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
//  x_ = x_ + K * (z - z_pred);
  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
