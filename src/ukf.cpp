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
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.0;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 1.0;

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
    // on incarnation, UKF is not initialized
    is_initialized_ = false;

    // start at 0 time step
    time_us_ = 0.;

    // set the state dimension
    n_x_ = 5;

    // update the augmented state dimension by ny_k and ny_phi_dot_dot
    n_aug_ = n_x_ + 2;

    // calculate size of the sigma point prediction matrix
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    // calculate size of the weights vector based on the sizes
    weights_ = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Make sure you switch between lidar and radar
    * measurements.
    */
    if(is_initialized_ == false) {

        switch(meas_package.sensor_type_) {
            case MeasurementPackage::RADAR:
            {

                // read values
                double rho = meas_package.raw_measurements_(0);
                double phi = meas_package.raw_measurements_(1);
                double rho_dot = meas_package.raw_measurements_(2);

                // calculate velocity
                double vx = rho_dot * cos(phi);
                double vy = rho_dot * sin(phi);
                double v = sqrt(vx * vx + vy * vy);

                // set the state vector
                x_ <<   std::min(rho * cos(phi), 0.0001),   // px
                        std::min(rho * sin(phi), 0.0001),   // üy
                        v,
                        rho,
                        rho_dot;

                // init covariance matrix with laser data
                P_ <<   std::pow(std_radr_, 2),     0,      0,      0,      0,
                        0,      std::pow(std_radr_, 2),     0,      0,      0,
                        0,      0,      std::pow(std_radrd_, 2),    0,      0,
                        0,      0,      0,      std_radphi_,        0,
                        0,      0,      0,      0,      std_radphi_;

                break;
            }

            case MeasurementPackage::LASER:
            {

                // can set values for position directly
                x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0., 0., 0.;

                // init covariance matrix
                P_ <<   std::pow(std_laspx_, 2),    0,      0,      0,      0,
                        0,      std::pow(std_laspy_, 2),    0,      0,      0,
                        0,      0,      1,      0,      0,
                        0,      0,      0,      1,      0,
                        0,      0,      0,      0,      1;
                break;
            }

            default:
                std::cout << "Error - unknown sensor type " << meas_package.sensor_type_ << std::endl;
                return;
        }

        // save measurement timestamp
        time_us_ = meas_package.timestamp_;

        // mark UKF as initialized
        is_initialized_ = true;

        // no predition as freshly initialized
        return;
    }

    // for prediction step, the time difference is needed - in seconds!
    double dT = (meas_package.timestamp_ - time_us_) / 1000000.0;

    // save measurement timestamp
    time_us_ = meas_package.timestamp_;

    // Prediction step with newly data
    Prediction(dT);

    // Update step -- need to take care what the measurement's source is
    switch(meas_package.sensor_type_) {
        case MeasurementPackage::RADAR:
            UpdateRadar(meas_package);
            break;

        case MeasurementPackage::LASER:
            UpdateLidar(meas_package);
            break;

        default:
            std::cout << "Error - unknown sensor type " << meas_package.sensor_type_ << std::endl;
            return;
    }
}

void UKF::Prediction(double delta_t) {
    /**
    * TODO: Complete this function! Estimate the object's location.
    * Modify the state vector, x_. Predict sigma points, the state,
    * and the state covariance matrix.
    */


    // generate sigma points over n_x_ dimension
    lambda_ = 3 - n_x_;

    // sigma point matrix creation
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

    // calculate square root of P
    MatrixXd A = P_.llt().matrixL();

    // set first coloumn of sigma point matrix with state vector
    Xsig.col(0) = x_;

    // fill remaining coloumns with calculated values
    for(int i = 0; i < n_x_; ++i) {
        Xsig.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
        Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
    }

    // spreading for augmentation
    lambda_ = 3 - n_aug_;

    // create augmented mean vector, augmented covariance matrix and sigma point matrix
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // create the augmented mean state
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0;
    x_aug(n_x_ + 1) = 0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std::pow(std_a_, 2);
    P_aug(n_x_ + 1, n_x_ + 1) = std::pow(std_yawdd_, 2);

    // create square root matrix again
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points - first col is state vector
    Xsig_aug.col(0) = x_aug;

    // calculate rest
    for(int i = 0; i < n_aug_; ++i) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    // start sigma point prediction
    VectorXd old = VectorXd(n_x_);
    VectorXd pred = VectorXd(n_x_);
    VectorXd noise = VectorXd(n_x_);

    for(int i = 0; i < 2 * n_aug_ + 1; i++) {

        VectorXd currCol    = Xsig_aug.col(i);
        double px           = currCol(0);
        double py           = currCol(1);
        double v            = currCol(2);
        double yaw          = currCol(3);
        double yawd         = currCol(4);
        double nu_a         = currCol(5);
        double nu_yawdd     = currCol(6);

        VectorXd old = currCol.head(5);

        // predict new position and velocity based on motion model
        if(yawd > .001) {
            double div = v / yawd;
            double vel = yaw + yawd * delta_t;

            pred << div * (sin(vel) - sin(yaw)),
                    div * (cos(yaw) - cos(vel)),
                   0,
                   yawd * delta_t,
                   0;
        } else {
            pred << v * cos(yaw) * delta_t,
                    v * sin(yaw) * delta_t,
                    0,
                    yawd * delta_t,
                    0;
        }

        // add noise to the prediction
        double dT_p2 = 0.5 * std::pow(delta_t, 2);

        noise <<    dT_p2 * cos(yaw) * nu_a,
                    dT_p2 * sin(yaw) * nu_a,
                    delta_t * nu_a,
                    dT_p2 * nu_yawdd,
                    delta_t * nu_yawdd;

        // add predicted sigma points
        Xsig_pred_.col(i) << old + pred + noise;
    }

    // prepare prediction state vector and covariance matrix
    VectorXd x_pred = VectorXd(n_x_);
    x_pred.fill(0.);

    MatrixXd P_pred = MatrixXd(n_x_, n_x_);
    P_pred.fill(0.);

    // set weights for sate mean calculation
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        // calculate weights
        if(i == 0)
            weights_(0) = lambda_ / (lambda_ + n_aug_);
        else
            weights_(i) = 0.5 / (n_aug_ + lambda_);

        // predict state mean
        x_pred += weights_(i) * Xsig_pred_.col(i);
    }

    // predict covariance matrix
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_pred;

        // angle normalization
        while(x_diff(3) > M_PI)
            x_diff(3) -= 2. * M_PI;

        while(x_diff(3) < -M_PI)
            x_diff(3) += 2. * M_PI;

        P_pred += weights_(i) * x_diff * x_diff.transpose();
    }

    x_ = x_pred;
    P_ = P_pred;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use lidar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the lidar NIS, if desired.
    */

    // set measurement state dimension --> lidar has px and py
    int n_z = 2;

    // create sigma point matrix, mean predicted measurement and covariance matrix in
    // measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.);

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.);

    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.);

    // transform sigma points into measurement space
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        // values for further calculation
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);

        // measurement model
        Zsig(0, i) = px;
        Zsig(1, i) = py;

        // mean predicted measurement
        z_pred += weights_(i) * Zsig.col(i);
    }

    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise to covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R <<    std::pow(std_laspx_, 2),    0,
            0,      std::pow(std_laspy_, 2);

    S += R;

    // create vector for incomming measurement
    VectorXd z = VectorXd(n_z);

    z <<    meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1);

    // create cross correlation matrix Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate the cross correlation matrix
    Tc.fill(0.);

    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        VectorXd z_diff = Zsig.col(i) - z_pred;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        while(x_diff(3) > M_PI)
            x_diff(3) -= 2. * M_PI;

        while(x_diff(3) < -M_PI)
            x_diff(3) += 2. * M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate the Kalman gain
    MatrixXd K = Tc * S.inverse();

    // residual
    VectorXd z_diff = z - z_pred;

    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use radar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the radar NIS, if desired.
    */

    // set measurement state dimension --> radar has r, phi, r_dot
    int n_z = 3;

    // create sigma point matrix, mean predicted measurement and covariance matrix in
    // measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.);

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.);

    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.);

    // transform sigma points into measurement space
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        // values for further calculation
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double vx = v * cos(yaw);
        double vy = v * sin(yaw);

        // measurement model
        Zsig(0, i) = sqrt(std::pow(p_x, 2) + std::pow(p_y, 2));     // r
        Zsig(1, i) = atan2(p_y, p_x);                               // phi
        Zsig(2, i) = (p_x * vx + p_y * vy) / Zsig(0, i);            // r_dot

        // mean predicted measurement
        z_pred += weights_(i) * Zsig.col(i);
    }

    //calculate measurement covariance matrix S
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        // calculate residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        while (z_diff(1) > M_PI)
            z_diff(1) -= 2. * M_PI;

        while (z_diff(1) < - M_PI)
            z_diff(1) += 2. * M_PI;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise to covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R <<    std::pow(std_radr_, 2),     0,      0,
            0,      std::pow(std_radphi_, 2),   0,
            0,      0,      std::pow(std_radrd_, 2);

    S += R;

    // create vector for incomming measurement
    VectorXd z = VectorXd(n_z);

    z <<    meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            meas_package.raw_measurements_(2);

    // create cross correlation matrix Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    // calculate the cross correlation matrix
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        VectorXd z_diff = Zsig.col(i) - z_pred;

        //normalize angles
        while (z_diff(1) > M_PI)
            z_diff(1) -= 2. * M_PI;

        while (z_diff(1) < -M_PI)
            z_diff(1) += 2. * M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        while(x_diff(3) > M_PI)
            x_diff(3) -= 2. * M_PI;

        while(x_diff(3) < -M_PI)
            x_diff(3) += 2. * M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate the Kalman gain
    MatrixXd K = Tc * S.inverse();

    // residual
    VectorXd z_diff = z - z_pred;

    // angle normalization
    while(z_diff(1) > M_PI)
        z_diff(1) -= 2. * M_PI;

    while(z_diff(1) < -M_PI)
        z_diff(1) += 2. * M_PI;


    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}
