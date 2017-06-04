#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define _MY_DEBUG
#undef _MY_DEBUG

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // TODO: Come back to this later

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // TODO: Come back to this later

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* State dimension
  time_us_ = 0;
  is_initialized_ = false;
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  lambda_ = 3 -n_aug_;
  NIS_laser_ = NIS_radar_ = 0;

  P_ << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1;

  int num_sig_pts = 2 * n_aug_ + 1;
  weights_ = VectorXd(num_sig_pts);
  for (int i=0; i < num_sig_pts; i++){
  		if (i == 0){
  			weights_(i) = lambda_/(lambda_+n_aug_);
  		}else{
  			weights_(i) = 1/(2*(lambda_+n_aug_));
  		}
  	}

}

UKF::~UKF() {}

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

	if (!is_initialized_) {
		// first measurement
		this->x_ = VectorXd(5);
		this->x_ << 1, 1, 1, 1, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		  /**
		  Convert radar from polar to cartesian coordinates and initialize state.
		  */
			float x = meas_package.raw_measurements_[0] * cos((double)meas_package.raw_measurements_[1]);
			float y = meas_package.raw_measurements_[0] * sin((double)meas_package.raw_measurements_[1]);
			this->x_ << x, y, 0, 0, 0;
			time_us_ = meas_package.timestamp_;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		  /**
		  Initialize state.
		  */
			this->x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
			time_us_ = meas_package.timestamp_;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	  }
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;


	this->Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		this->UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		this->UpdateLidar(meas_package);
	}

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

	int num_sig_pts = 2 * n_aug_ + 1;

	// Create Augmented mean vector
	VectorXd x_aug  = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	MatrixXd Q = MatrixXd(2,2);
	Q << this->std_a_ * this->std_a_, 0, 0, this->std_yawdd_*this->std_yawdd_;

	//
	x_aug.head(this->n_x_) = this->x_;
	x_aug(n_x_+0) = 0;
	x_aug(n_x_+1) = 0;


	P_aug = MatrixXd::Zero(n_aug_, n_aug_);

	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug.bottomRightCorner(2,2) = Q;


	//calculate square root of P
	MatrixXd A = P_aug.llt().matrixL();



	//set first column of sigma point matrix
	Xsig_aug.col(0)  = x_aug;

	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
	}


#ifdef _MY_DEBUG
	//print result
	std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
#endif

	// Now predict the sigma points
	this->Xsig_pred_= MatrixXd(n_x_, num_sig_pts);

	for(int i=0; i < num_sig_pts; i++){

		double v = Xsig_aug(2, i);
		double psi = Xsig_aug(3, i);
		double psi_dot = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_psi = Xsig_aug(6, i);
		//std::cerr << i << " : " << nu_a << ":"<<(1.0/2.0)*(delta_t * delta_t) << endl; // << " : " <<  cos(psi) << " : " << nu_a << " : " << ((1/2)*(delta_t * delta_t)*cos(psi)*nu_a) << std::endl;

		if (psi_dot != 0){
			Xsig_pred_(0, i) = Xsig_aug(0, i) +  (v/psi_dot)*( sin(psi + psi_dot*delta_t) - sin(psi)) + (1.0/2.0)*(delta_t * delta_t)*cos(psi)*nu_a;
			Xsig_pred_(1, i) = Xsig_aug(1, i) +  (v/psi_dot)*( -cos(psi + psi_dot*delta_t) + cos(psi)) + (1.0/2.0)*(delta_t * delta_t)*sin(psi)*nu_a;
			Xsig_pred_(2, i) = Xsig_aug(2, i) + 0 + delta_t * nu_a;
			Xsig_pred_(3, i) = Xsig_aug(3, i) + psi_dot * delta_t + (1.0/2.0)*(delta_t*delta_t)*nu_psi;
			Xsig_pred_(4, i) = Xsig_aug(4, i) + 0 + delta_t * nu_psi;
		}else{
			Xsig_pred_(0, i) = Xsig_aug(0, i) +  v*cos(psi)*delta_t +  (1.0/2.0)*(delta_t * delta_t)*cos(psi)*nu_a;
			Xsig_pred_(1, i) = Xsig_aug(1, i) +  v*sin(psi)*delta_t +  (1.0/2.0)*(delta_t * delta_t)*sin(psi)*nu_a;;
			Xsig_pred_(2, i) = Xsig_aug(2, i) + 0 + delta_t * nu_a;
			Xsig_pred_(3, i) = Xsig_aug(3, i) + psi_dot * delta_t + (1.0/2.0)*(delta_t*delta_t)*nu_psi;
			Xsig_pred_(4, i) = Xsig_aug(4, i) + 0 + delta_t * nu_psi;
		}

	}
#ifdef _MY_DEBUG
	//print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
#endif

	// Now predict mean and covariance

	VectorXd x = VectorXd::Zero(n_x_);
	MatrixXd P = MatrixXd::Zero(n_x_, n_x_);


	for (int i=0; i < num_sig_pts; i++){
		x += weights_(i) * Xsig_pred_.col(i);
	}

	for (int i=0; i < num_sig_pts; i++){
		VectorXd diff = Xsig_pred_.col(i) - x;
		P += weights_(i) * diff * diff.transpose();
	}
#ifdef _MY_DEBUG
	//print result
	std::cout << "Predicted state" << std::endl;
	std::cout << x << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P << std::endl;
#endif

	// Now overwrite the results and thus prediction is completed
	this->x_ = x;
	this->P_ = P;



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
// First change prediction to measurement space

	int n_z = 2; // Number of radar measurement dimensions

	double lambda = 3 - n_aug_;
	int num_sig_pts = 2*n_aug_+1;

	// already initialized
	std_laspx_ = 0.15;
	std_laspy_ = 0.15;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, num_sig_pts);
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z, n_z);

	// First transform to measurement space
	for (int i=0; i<num_sig_pts; i++){
		Zsig(0, i) = Xsig_pred_(0, i);
		Zsig(1, i) = Xsig_pred_(1, i);
	}

	// Now calculate mean
	z_pred.fill(0);
	for (int i=0; i<num_sig_pts; i++){
		z_pred += weights_(i) * Zsig.col(i);
	}

	// And variance
	MatrixXd R = MatrixXd::Zero(2, 2);
	R(0, 0) = std_laspx_ * std_laspx_;
	R(1, 1) = std_laspy_ * std_laspy_;

	S.fill(0);
	for (int i=0; i<num_sig_pts; i++){
		VectorXd diff = Zsig.col(i) - z_pred;
		S += weights_(i) * diff * diff.transpose();
	}
	S += R;
#ifdef _MY_DEBUG
	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;
#endif

	//
	// Now update
	//

  // New measurement at k+1
	VectorXd z = VectorXd(n_z);
	z <<
	  meas_package.raw_measurements_[0],   //px
	  meas_package.raw_measurements_[1];   //py


	// Correlation Tc of predicted space and measurement space
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	// Now do the work
	// 1. Calculate correlation matrix Tc
	for (int i=0; i < num_sig_pts; i++){
	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	  VectorXd z_diff = Zsig.col(i) - z_pred;
	  Tc += weights_(i) * (x_diff) * (z_diff).transpose();
	}
	// 2. Calculate Kalman gain
	MatrixXd K = Tc*S.inverse();

	// 3. Update State
	x_ = x_ + K*(z-z_pred);

	// 4. Update covariance
	P_ = P_ - K*S*K.transpose();

#ifdef _MY_DEBUG
	//print result
	std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
#endif
	VectorXd Zdiff = z-z_pred;
	NIS_laser_ = (Zdiff).transpose()*S.inverse()*(Zdiff);

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

	// First change prediction to measurement space

	int n_z = 3; // Number of radar measurement dimensions

	double lambda = 3 - n_aug_;
	int num_sig_pts = 2*n_aug_+1;


	double std_rad_rho = std_radr_;
	double std_rad_phi = std_radphi_;
	double std_rad_rho_dot = std_radrd_;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, num_sig_pts);

	VectorXd z_pred = VectorXd(n_z);

	MatrixXd S = MatrixXd(n_z, n_z);

	// First transform to measurement space
	for (int i=0; i<num_sig_pts; i++){
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);
		//double psi_dot = Xsig_pred(4, i); // not needed

		double rho = sqrt(px*px + py*py);
		double phi = atan2(py, px);
		while (phi > M_PI) phi-=2.*M_PI;
		while (phi < -M_PI) phi+=2.*M_PI;
		double rho_dot = (px * cos(psi)* v + py * sin(psi) * v)/rho;
		Zsig(0, i) = rho;
		Zsig(1, i) = phi;
		Zsig(2, i) = rho_dot;

	}

	// Now calculate mean
	z_pred.fill(0);
	for (int i=0; i<num_sig_pts; i++){
		z_pred += weights_(i) * Zsig.col(i);
	}

	// And variance
	MatrixXd R = MatrixXd::Zero(3, 3);
	R(0, 0) = std_rad_rho * std_rad_rho;
	R(1, 1) = std_rad_phi * std_rad_phi;
	R(2, 2) = std_rad_rho_dot * std_rad_rho_dot;

	S.fill(0);
	for (int i=0; i<num_sig_pts; i++){
		VectorXd diff = Zsig.col(i) - z_pred;
		while (diff(1)> M_PI) diff(1)-=2.*M_PI;
		while (diff(1)<-M_PI) diff(1)+=2.*M_PI;
		S += weights_(i) * diff * diff.transpose();
	}
	S += R;
#ifdef _MY_DEBUG
	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;
#endif
	//
	// Now update
	//

		  // New measurement at k+1
	  VectorXd z = VectorXd(n_z);
	  z <<
		  meas_package.raw_measurements_[0],   //rho in m
		  meas_package.raw_measurements_[1],   //phi in rad
		  meas_package.raw_measurements_[2];   //rho_dot in m/s

	  // Correlation Tc of predicted space and measurement space
	  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	  // Now do the work
	  // 1. Calculate correlation matrix Tc
	  for (int i=0; i < num_sig_pts; i++){
          VectorXd x_diff = Xsig_pred_.col(i) - x_;
		  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		  VectorXd z_diff = Zsig.col(i) - z_pred;
		  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		  Tc += weights_(i) * (x_diff) * (z_diff).transpose();
	  }

	  // 2. Calculate Kalman gain
	  MatrixXd K = Tc*S.inverse();

	  // 3. Update State
	  x_ = x_ + K*(z-z_pred);

	  // 4. Update covariance
	  P_ = P_ - K*S*K.transpose();

	  VectorXd Zdiff = z-z_pred;
	  NIS_radar_ = (Zdiff).transpose()*S.inverse()*(Zdiff);
#ifdef _MY_DEBUG
	  //print result
	  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
#endif

}
