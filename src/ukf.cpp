#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  n_x = 5;
  n_aug = 7;
  
  VectorXd weights = VectorXd(2*n_aug+1);
  
  
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
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
    cout << "UKF: " << endl;
    x_ << 0,0,1,1,0.1;
	
    P_ << 0.1,0,0,0,0,
		0,0.1,0,0,0,
		0,0,1,0,0,
		0,0,0,1,0,
		0,0,0,0,1;
	
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      
		double rho = meas_package.raw_measurements_[0];
		double fi = meas_package.raw_measurements_[1];
		double yaw = meas_package.raw_measurements_[2];
		x_(0) = rho*cos(fi);
		x_(1) = rho*sin(fi);	 
		}else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {      
		x_(0) = meas_package.raw_measurements_[0];
		x_(1) = meas_package.raw_measurements_[1];
		}
		

    // done initializing, no need to predict or update
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;	
    return;
	}
		
	double delta = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;	
	
	Prediction(delta);
	
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
	UpdateRadar(meas_package);
		}
	
		
	if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
		UpdateLidar(meas_package);
		}
	cout<< "x_    "<<x_<<endl;
	cout<< "P_     "<<P_<<endl;
	
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
  double dt = delta_t;
  //define spreading parameter
  lambda = 3 - n_x;
  lambda_aug = 3 - n_aug;
	
	MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);
	
	MatrixXd sqrtP = P_.llt().matrixL();
	
    Xsig.col(0) = x_;
	
	for (int i = 0; i < n_x; i++){
        Xsig.col(i+1) = x_ + sqrt(lambda + n_x) * sqrtP.col(i);
        Xsig.col(i+1+n_x) = x_ - sqrt(lambda + n_x) * sqrtP.col(i);
	}
	
//################# Augmented Sigma points

	//define spreading parameter
  
		
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd sqrtP_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  
  
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_aug + n_aug) * sqrtP_aug.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_aug + n_aug) * sqrtP_aug.col(i);
  }
  
// ########### Augmented to Sigma transformation

  //create matrix with predicted sigma points as columns
  Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
//predict sigma points
  
  for (int i=0; i<2*n_aug+1; i++){
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double w = Xsig_aug(3,i);
      double wd = Xsig_aug(4,i);
      double nuv = Xsig_aug(5,i);
      double nuw = Xsig_aug(6,i);
      
      
      if (fabs(wd) > 0.001){
      px +=  v*(sin(w+wd*dt) - sin(w))/wd + 0.5*pow(dt,2)*cos(w)*nuv;
      py +=  v*(-cos(w+wd*dt) + cos(w))/wd + 0.5*pow(dt,2)*sin(w)*nuv;      
      }else{
      px +=  v*cos(w)*dt + 0.5*pow(dt,2)*cos(w)*nuv;
      py +=  v*sin(w)*dt + 0.5*pow(dt,2)*sin(w)*nuv;
      }
	  
	  v += dt*nuv;
      w += wd*dt + 0.5*pow(dt,2)*nuw;
      wd += dt*nuw;
      	     
    Xsig_pred(0,i) = px;
    Xsig_pred(1,i) = py;
    Xsig_pred(2,i) = v;
    Xsig_pred(3,i) = w;
    Xsig_pred(4,i) = wd;
	}
	
	//create vector for weights
  weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda_aug/(lambda_aug+n_aug);
 
  //##### weight calculation  
  weights(0) = weight_0;
  
  for (int i = 1 ; i<2*n_aug+1 ; i++){
	  double weight_ = 0.5/(lambda_aug+n_aug);
      weights(i)  = weight_;
  }
  
  x_.fill(0.0);
  
  for(int i = 0; i<2*n_aug+1; i++){
      x_ +=  weights(i) * Xsig_pred.col(i);
  }
  P_.fill(0.0);
  for(int i = 0; i<2*n_aug+1; i++){
      VectorXd lk = Xsig_pred.col(i) - x_;
	  
	  if (lk(3)> M_PI)
		lk(3) -= 2.*M_PI;
	
	  if (lk(3)<  -M_PI)
		lk(3) += 2.*M_PI;
	  
      P_ += weights(i) * lk * lk.transpose();
  }
  
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
  cout<<"Lidar"<<endl;
    int n_l = 2;
	VectorXd z = VectorXd(n_l);
	z << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1];
	
	
	MatrixXd Zsig = MatrixXd(n_l, 2 * n_aug + 1);
	VectorXd z_pred = VectorXd(n_l);
	MatrixXd S = MatrixXd(n_l,n_l);
	
     for (int i = 0; i<2 * n_aug + 1; i++){
      Zsig(0,i) = Xsig_pred(0,i);
      Zsig(1,i) = Xsig_pred(1,i);
		}
  //calculate mean predicted measurement
	z_pred.fill(0.0);
	 for (int i = 0; i<2 * n_aug + 1; i++){
        z_pred += weights(i)*Zsig.col(i);
		}
		cout<<"beafore S"<<endl;
		S.fill(0.0);
  //calculate innovation covariance matrix S
	 for (int i = 0; i<2 * n_aug + 1; i++){
        VectorXd kl = Zsig.col(i) - z_pred;
        S += weights(i) * kl * kl.transpose();
		}
		
	MatrixXd R = MatrixXd(n_l,n_l);
             R << std_laspx_*std_laspx_ , 0 ,
				  0 ,std_laspy_*std_laspy_;
	S += R;
  cout<<"after S"<<endl;
	MatrixXd Tc = MatrixXd(n_x, n_l);

	Tc.fill(0.0);
  //calculate cross correlation matrix
	 for(int i = 0; i< 2 * n_aug + 1 ; i++){
		VectorXd Xk = Xsig_pred.col(i) - x_;
		/*if (Xk(3)> M_PI)
		Xk(3) -= 2.*M_PI;
	
		if (Xk(3) < -M_PI)
		Xk(3) += 2.*M_PI;
		*/
		VectorXd Zk = Zsig.col(i) - z_pred;
		Tc += weights(i) * Xk * Zk.transpose();
		}
  
  //calculate Kalman gain K;
  cout<<"b4 kalman"<<endl;
	MatrixXd Kalman = Tc * S.inverse();
  //update state mean and covariance matrix
    x_ += Kalman * (z - z_pred);
    P_ -= Kalman * S * Kalman.transpose();  
	cout<<"sdfgsf"<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_. **/
  cout<<"radar"<<endl;
  
  int n_z = 3;
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];
	
	MatrixXd Zsigr = MatrixXd(n_z, 2 * n_aug + 1);
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z,n_z);
	
     for (int i = 0; i < 2 * n_aug + 1; i++){
		 
		double p_x = Xsig_pred(0,i);
		double p_y = Xsig_pred(1,i);
		double v_  = Xsig_pred(2,i);
		double yaw_ = Xsig_pred(3,i);

		double v1 = cos(yaw_)*v_;
		double v2 = sin(yaw_)*v_;
		
		/*if (p_x<0.001 ){
			p_x = 0.001;
		}
		
		if (p_y<0.001){
		p_y = 0.001;
		}*/
		
	    Zsigr(0,i) = sqrt(p_x * p_x + p_y * p_y);                        //r
		Zsigr(1,i) = atan2(p_y,p_x);
		Zsigr(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}
	
  //calculate mean predicted measurement
		
		 z_pred.fill(0.0);
	 for (int i = 0; i<2 * n_aug + 1; i++){
        z_pred += weights(i)*Zsigr.col(i);
		}
  
		S.fill(0.0);
  //calculate innovation covariance matrix S
	 for (int i = 0; i<2 * n_aug + 1; i++){
        VectorXd kl = Zsigr.col(i) - z_pred;
		
		if (kl(1)> M_PI)
			kl(1)-=2.*M_PI;
		if (kl(1)<-M_PI)
			kl(1)+=2.*M_PI;
		
        S += weights(i)*kl*kl.transpose();
		}
		
	MatrixXd R = MatrixXd(n_z,n_z);
             R << std_radr_*std_radr_ , 0 , 0,
				  0 ,std_radphi_*std_radphi_ ,0,
				  0,0, std_radrd_*std_radrd_;
	S += R;
  
	MatrixXd Tc = MatrixXd(n_x, n_z);

	Tc.fill(0.0);
  //calculate cross correlation matrix
	 for(int i = 0; i< 2 * n_aug + 1 ; i++){
		VectorXd Xk = Xsig_pred.col(i) - x_;
		if (Xk(3)> M_PI)
			Xk(3)-=2.*M_PI;
		if (Xk(3)<-M_PI)
			Xk(3)+=2.*M_PI;
		
		VectorXd Zk = Zsigr.col(i) - z_pred;
		if (Zk(1)> M_PI)
			Zk(1)-=2.*M_PI;
		if (Zk(1)<-M_PI)
			Zk(1)+=2.*M_PI;
		
		Tc += weights(i) * Xk * Zk.transpose();
		}
  
  //calculate Kalman gain K;
	VectorXd z_diff = z - z_pred;

  //angle normalization
  if (z_diff(1)> M_PI) 
	  z_diff(1)-=2.*M_PI;
  if (z_diff(1)<-M_PI) 
	  z_diff(1)+=2.*M_PI;
  
	MatrixXd Kalman = Tc * S.inverse();
  //update state mean and covariance matrix
    x_ += Kalman * z_diff;
    P_ -= Kalman * S * Kalman.transpose();   
  /*
  You'll also need to calculate the radar NIS.
  */
}

