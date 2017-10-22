#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  Init();
}

UKF::~UKF() {}

/**
 * initializtion function
 * 
 */
void UKF::Init()
{
   // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

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

   // dimension of state vector
   n_x_ = 5;
   // dimension of augmentation state vector
   n_aug_ = 7;
   // lambda
   lambda_ = 3 - n_x_;

   // initial state vector
   x_ = VectorXd(n_x_);
   // initial covariance matrix
   P_ = MatrixXd(n_x_, n_x_);
   P_ <<     1.0,   0,    0,   0,   0,
            0,    1.0,    0,    0,    0,
            0,    0,    1.0,    0,    0,
            0,    0,    0,    1.0,    0,
            0,    0,    0,    0,    1.0;
   // augmentation state vector
   x_aug_ = VectorXd(n_aug_);
   // augmentation covariance matrix
   P_aug_ = MatrixXd(n_aug_, n_aug_);

   weights_ = VectorXd(2 * n_aug_ + 1);
   // sigma points matrix
   Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);
   // augmentation sigma points matrix
   Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
   // predicted sigma points matrix
   Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

   //set weights
   lambda_ = 3 - n_aug_;
   weights_.fill(0.5/(lambda_+n_aug_));
   weights_(0)=lambda_/(lambda_+n_aug_);
   
   //set measurement dimension, radar can measure r, phi, and r_dot
   n_z_ = 3;
   //mean predicted measurement
   z_pred_ = VectorXd(n_z_);
   //create matrix for sigma points in measurement space
   Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
   //measurement covariance matrix S
   S_ = MatrixXd(n_z_,n_z_);

   //whether initialized x_
   is_initialized_=false;
   //timestamp
   previous_timestamp_=0;

}

void UKF::GenerateSigmaPoints() {

  lambda_= 3 - n_x_;

  //calculate square root of P_
  MatrixXd A_ = P_.llt().matrixL();

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig
  Xsig_.col(0)=x_;
  for(int i=0;i<n_x_;i++)
  {   
      Xsig_.col(i+1) = x_ + sqrt(lambda_ + n_x_)*A_.col(i);
      Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_)*A_.col(i);
  }

}

void UKF::AugmentedSigmaPoints() {


  lambda_ = 3 - n_aug_;
 
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
}

void UKF::SigmaPointPrediction(double delta_t) {


//predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}

void UKF::PredictMeanAndCovariance() {


  lambda_ = 3 - n_aug_;

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);

  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}


void UKF::PredictRadarMeasurement() {


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd(n_z_,n_z_);
  R_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_ = S_ + R_;

 
}

void UKF::PredictLaserMeasurement() {


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig_(0,i) = p_x; //px
    Zsig_(1,i) = p_y; //py
    Zsig_(2,i) = 0;   //0
  }

  //mean predicted measurement
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd(n_z_,n_z_);
  R_ <<    std_laspx_*std_laspx_, 0, 0,
          0, std_laspy_*std_laspy_, 0,
          0, 0, 0.001;
  S_ = S_ + R_;

}


void UKF::UpdateState(VectorXd z) {


  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K_ = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff;
  P_ = P_ - K_*S_*K_.transpose();

  //cout << "z = " << z << endl;
  //cout << "Tc_ = " << Tc_ << endl;
  //cout << "S_.inverse() = " << S_.inverse() << endl;
  //cout << "K_ = " << K_ << endl;
  //cout << "z_diff = " << z_diff << endl;
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}


double UKF::CalcuNIS(VectorXd z) {

  //residual
  VectorXd z_diff = z - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //calculate the epsilong(NIS)
  double epsilong = z_diff.transpose()*S_.inverse()*z_diff;
  return epsilong;
}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
    */
    // first measurement
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)&&(use_radar_ == true)) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_(0)=meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      x_(1)=meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
      x_(2)=meas_package.raw_measurements_[2];
      x_(3)=meas_package.raw_measurements_[1];
      x_(4)=0;
     
    }
    else if ((meas_package.sensor_type_ == MeasurementPackage::LASER)&&(use_laser_ == true)) {
      /**
      Initialize state.
      */
      x_(0)=meas_package.raw_measurements_[0];
      x_(1)=meas_package.raw_measurements_[1];
      x_(2)=0;
      x_(3)=0;
      x_(4)=0;
    }
    else
    {
      return;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)&&(use_radar_ == true)) {

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //dt - expressed in seconds
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

    // Radar updates   
    UpdateRadar(meas_package);

  } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER)&&(use_laser_ == true)){
   /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //dt - expressed in seconds
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
   // Laser updates

    UpdateLidar(meas_package);
  }
  else{
   return;
  }
  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //cout << "x_ = " << x_ << endl;
  GenerateSigmaPoints();
  //cout << "Xsig_ = " << Xsig_ << endl;
  AugmentedSigmaPoints();
  //cout << "Xsig_aug_ = " << Xsig_aug_ << endl;
  SigmaPointPrediction(delta_t);
  //cout << "Xsig_pred_ = " << Xsig_pred_ << endl;
  PredictMeanAndCovariance();
  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  PredictLaserMeasurement();

  VectorXd z=VectorXd(n_z_);
  z<<meas_package.raw_measurements_[0],
     meas_package.raw_measurements_[1],
     0;
  UpdateState(z);

  NIS_laser_ = CalcuNIS(z);
  cout << "NIS_laser_ = " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  PredictRadarMeasurement();
  
  VectorXd z=VectorXd(n_z_);
  z<<meas_package.raw_measurements_[0],
     meas_package.raw_measurements_[1],
     meas_package.raw_measurements_[2];
  UpdateState(z);

  NIS_radar_ = CalcuNIS(z);
  cout << "NIS_radar_ = " << NIS_radar_ << endl;
}
