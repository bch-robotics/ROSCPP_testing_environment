// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <memory>
#include <vector>
#include <bitset>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <unordered_map>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <ctime>
#include <algorithm>
#include <stdexcept>
#include <cfenv>
#include <chrono>
#include <array>

typedef Eigen::Matrix< long double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > Mat;

class MedianFilter {
    private:
        std::vector<double> window;
        size_t windowSize;
        size_t middleIndex;

    public:
        MedianFilter(size_t size) : windowSize(size), middleIndex(size / 2) {
            window.reserve(windowSize);
        }

        double calculate(double newSample) {
            if (window.size() < windowSize) {
                window.push_back(newSample);
            }
            else {
                window[windowSize - 1] = newSample;
            }

            std::vector<double> sortedWindow = window;
            std::nth_element(sortedWindow.begin(), sortedWindow.begin() + middleIndex, sortedWindow.end());

            return sortedWindow[middleIndex];
        }
};


class ExponentialMovingAverage
{
private:
    double alpha;
    double storedValue;
    bool uninitialized;

public:
    ExponentialMovingAverage(double alpha) : storedValue(0), uninitialized(true)
    {
        if (alpha <= 0 || alpha > 1)
        {
            throw std::invalid_argument("Alpha must be in the range (0, 1]");
        }
        this->alpha = alpha;
    }

    double calculate(double toAdd)
    {
        if (toAdd > std::numeric_limits<double>::max() || toAdd < std::numeric_limits<double>::lowest())
        {
            throw std::overflow_error("Numerical overflow or underflow in added value");
        }

        if (this->uninitialized)
        {
            this->storedValue = toAdd;
            this->uninitialized = false;
            std::cout << "hi" << std::endl;
        }

        std::cout << "hi" << this->storedValue << std::endl;
        std::cout << "hey" << toAdd << std::endl;
        std::cout << "yo" << this->alpha << std::endl;

        double ema = this->storedValue + this->alpha * (toAdd - this->storedValue);
        if (ema > std::numeric_limits<double>::max() || ema < std::numeric_limits<double>::lowest())
        {
            throw std::overflow_error("Numerical overflow or underflow in EMA calculation");
        }

        this->storedValue = ema;
        return ema;
    }

    double getVal()
    {
        return this->storedValue;
    }
};

class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;

public:
    void tic() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    double toc() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        return duration.count() / 1e6; // Convert to seconds
    }
};

int numTubes = 3;
int numConfigParams = 3;
int numConfigParams_ = 3;
int numTotalTubes = 3;
int numTendonsPerTube = 2;
int nOfTendonMotors = 6;
int nOfDrivenMotors = 9;

const double PI = 3.14159265358979323846;
const long double EPSILON = 1e-6;
const int MAX_ITERATIONS = 100;

Vec automated_velocity(3);
Eigen::Matrix<long double, 2, 2, Eigen::RowMajor> tendonGeometry[3]; // creates a 3x2x2 matrix ie 3 2x2 matrices
Eigen::Matrix<long double, 2, 2, Eigen::RowMajor> invTendonGeometry[3]; // creates a 3x2x2 matrix ie 3 2x2 matrices

Vec theta_init_(numConfigParams); // initial bending angles of each sheath
Vec phi_init_(numConfigParams); // initial bending plane offsets of each sheath
Vec length_init_(numConfigParams); // initial length of the steerable section of each sheath in mm
Vec delta_(numConfigParams); // radius (in meters) of cross sectional circle centered at sheath centerline and passes through the tendon locations
long double lengthLowLimit = 1.0;  // in mm
long double lengthHighLimit = 50.0; // in mm
long double minAllowedRadii = 12.0; // in mm, calculated using length = radius * theta (in radians) where theta is 90 degrees
Vec maxJointSpeed(9);
Vec jointMidPoints(9);
int outOfWorkspace = 0;
long double previous_time_ = 0.0;
long double deltaT = 0.1;
long double redundancyResolutionOn = 1.0;
bool useOnlyPositionJacobian = true;
bool applyJointLimits = true;
/*  operating modes
    0: only position (3 DOF) with 5 DOF (no phi1 and sheath3),
    1: position and orientation with no roll (5 DOF) with 6 DOF robot (no phi1, theta3, and phi3)
    2: position and orientation (6 DOF) with 7 DOF robot (no phi1 and theta3).
*/
int operatingMode = 1;

bool exp_start_time_saved_ = false;
long double exp_start_time_ = 0.0;
int unemployed = 0;

// FOR DEBUGGING AND SIM 
Vec q_sim(9);
Vec qDot_sim(9);
Vec pTip_init_sim;
Vec pTip_sim;
Vec displacement_sim(9);

template <typename T> long double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Mat genIdent3() {
    Mat toReturn(3, 3);
    toReturn << long double(1.0), long double(0.0), long double(0.0), long double(0.0), long double(1.0), long double(0.0), long double(0.0), long double(0.0), long double(1.0);
    return toReturn;
}

Mat genZero3() {
    Mat toReturn(3, 3);
    toReturn << long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0);
    return toReturn;
}

Vec genZeros(int size) {
    Vec toReturn(size);
    for (int i = 0; i < size; ++i) {
        toReturn(i) = long double(0.0);
    }
    return toReturn;
}

Mat genIdentN(int size) {
    Mat toReturn(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                toReturn(i, j) = long double(1.0);
            }
            else {
                toReturn(i, j) = long double(0.0);
            }
        }
    }
    return toReturn;
}

Mat genZeroN(int size) {
    Mat toReturn(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            toReturn(i, j) = long double(0.0);
        }
    }
    return toReturn;
}

Mat skew(Vec w)
{
    if (3 != w.size())
    {
        throw std::invalid_argument("Vector must be 3x1");
    }

    Mat R(3, 3);
    R = genZero3();

    R(0, 1) = -1.0 * w(2);
    R(0, 2) = w(1);
    R(1, 2) = -1.0 * w(0);

    R(1, 0) = w(2);
    R(2, 0) = -1.0 * w(1);
    R(2, 1) = w(0);

    return R;
}

Mat pinv(Mat A, double tol)
{
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //Eigen::VectorXd s = svd.singularValues();

    ////std::cout << "single ladies" << s << std::endl;

    //if (!s.isZero(0))
    //{
    //    // int r1 = (s.array() > tol).count() + 1; // since the last s is too small and this +1 allows it to exist, its inverse down the line creates big values. This is to ensure that r1 is at least 1 but it is now done as follows (but matlab's pinv if you look at it carefully does not do it. although there is +1):
    //    int r1 = (s.array() > tol).count(); // estimate effective rank
    //    r1 = std::max(r1, 1);               // Ensure that r1 is at least 1 // matlab does not ensure this

    //    Eigen::MatrixXd U = svd.matrixU();
    //    Eigen::MatrixXd V = svd.matrixV();

    //    // V.rightCols(V.cols() - r1).setZero();
    //    V.conservativeResize(V.rows(), r1);
    //    // U.rightCols(U.cols() - r1).setZero();
    //    U.conservativeResize(U.rows(), r1);
    //    s.conservativeResize(r1);

    //    s = s.array().inverse();
    //    return V * s.asDiagonal() * U.transpose();
    //}
    //else
    //{
    //    return Eigen::MatrixXd::Constant(A.cols(), A.rows(), std::numeric_limits<double>::quiet_NaN());
    //}

    auto toReturn = A.completeOrthogonalDecomposition();
    toReturn.setThreshold(tol); // sets the threshold for which values should be considered to be 0
    return Mat(toReturn.pseudoInverse());
}

long double deg2rad(long double deg) {
    return deg * M_PI / 180.0;
}

long double rad2deg(long double rad) {
    return rad * 180.0 / M_PI;
}

Mat getRotx(long double ang) {
    Mat toReturn(3, 3);
    toReturn << 1.0, 0.0, 0.0,
        0.0, std::cos(deg2rad(ang)), -1.0 * std::sin(deg2rad(ang)),
        0.0, std::sin(deg2rad(ang)), std::cos(deg2rad(ang));

    return toReturn;
}

Mat getRoty(long double ang) {
    Mat toReturn(3, 3);
    toReturn << std::cos(deg2rad(ang)), 0.0, std::sin(deg2rad(ang)),
        0.0, 1.0, 0.0,
        -1.0 * std::sin(deg2rad(ang)), 0.0, std::cos(deg2rad(ang));

    return toReturn;
}

Mat getRotz(long double ang) {
    Mat toReturn(3, 3);
    toReturn << std::cos(deg2rad(ang)), -1.0 * std::sin(deg2rad(ang)), 0.0,
        std::sin(deg2rad(ang)), std::cos(deg2rad(ang)), 0.0,
        0.0, 0.0, 1.0;

    return toReturn;
}

Mat getRotx_rad(long double ang_rad) {
    return getRotx(rad2deg(ang_rad));
}

Mat getRoty_rad(long double ang_rad) {
    return getRoty(rad2deg(ang_rad));
}

Mat getRotz_rad(long double ang_rad) {
    return getRotz(rad2deg(ang_rad));
}

std::vector<std::vector<std::string>> readCSV(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        data.push_back(row);
    }

    return data;
}

void readCSVDub(const std::string& filename, std::map<int, Eigen::Matrix3d>& toReturn) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            try {
                double value = std::stod(cell);
                row.push_back(value);
            }
            catch (const std::invalid_argument& e) {
                // Handle non-numeric data (e.g., headers)
                // For simplicity, we'll skip this cell
                continue;
            }
            catch (const std::out_of_range& e) {
                // Handle out of range errors
                std::cerr << "Error: Number out of range - " << cell << std::endl;
                continue;
            }
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    int othercount = 0;
    // Access data like a 2D vector
    for (const auto& row : data) {
        Eigen::Matrix3d input = Eigen::Matrix3d::Identity();
        int count = 0;
        for (const auto& value : row) {
            input(count) = value;
            count++;
        }
        toReturn[othercount] = input;
        othercount++;
    }
}

void writeCSV(const std::string& filename, const std::vector<Eigen::VectorXd>& data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}

void writeVectorToCSV(const std::string& filename, const Vec& vector) {
    std::ofstream file;

    //if (!file.is_open()) {
    //    throw std::runtime_error("Unable to open file: " + filename);
    //}

    file.open(filename, std::ios::app);

    file << vector.format(Eigen::IOFormat(
        Eigen::FullPrecision,
        Eigen::DontAlignCols,
        ",",
        "\n"
    ));
    file << "\n";
    file.close();
}

int find_element(std::vector<double> toSearch, double toFind) {
    int toReturn = -1;
    for (int i = 0; i < toSearch.size(); ++i) {
        toReturn = (toSearch[i] == toFind) ? i : toReturn;
    }
    return toReturn;
}

/*
    Calculates position and orientation of robot tip
    with respect to the global reference frame
*/
void forwardKinematics(Vec q,
    Vec& end_tip_pos,
    Mat& rotationMatrix) {

    int numSheaths = q.rows() / numConfigParams;

    Vec theta(numSheaths);
    Vec phi(numSheaths);
    Vec length(numSheaths);

    // store the config space params in separate vectors
    for (int i = 0; i < numSheaths; ++i) {
        theta(i) = (q(i * numConfigParams) == 0.0) ? std::numeric_limits<long double>::epsilon() : q(i * numConfigParams);
        phi(i) = q(i * numConfigParams + 1);
        length(i) = q(i * numConfigParams + 2);
    }

    // initialize rotation matrix and end tip position vector
    rotationMatrix = genIdent3();
    end_tip_pos = genZeros(3);

    for (int k = 0; k < numSheaths; ++k) {
        rotationMatrix = genIdent3();

        for (int j = 0; j < (k + 1); ++j) {

            if (j < k) { // calculate rotation matrix if not at final sheath
                rotationMatrix = rotationMatrix * (getRotz_rad(phi(j))) * (getRoty_rad(theta(j)));
            }
            else { // calculate end tip position of current sheath wrt previous frame
                Vec tip_pos_rel(3);
                tip_pos_rel << (length(j) / theta(j)) * (1.0 - std::cos(theta(j))),
                    0.0,
                    (length(j) / theta(j))* std::sin(theta(j));
                end_tip_pos = end_tip_pos + rotationMatrix * (getRotz_rad(phi(j))) * tip_pos_rel;
            }

            if (j == (numSheaths - 1)) { // calculate rotation matrix of final sheath wrt global frame
                rotationMatrix = rotationMatrix * (getRotz_rad(phi(j))) * (getRoty_rad(theta(j)));
            }
        }
    }
}

/*
    Differential kinematics calculation that returns
    a Jacobian J given configuration parameters q
*/
Mat generateJacobian(Vec q) {

    int numSheaths = q.rows() / numConfigParams;

    Vec theta(numSheaths);
    Vec phi(numSheaths);
    Vec length(numSheaths);

    // separate and store config params
    for (int i = 0; i < numSheaths; ++i) {
        theta(i) = q(i * numConfigParams);
        phi(i) = q(i * numConfigParams + 1);
        length(i) = q(i * numConfigParams + 2);
    }

    Mat J(6, numSheaths * numConfigParams); // Jacobian matrix, 6x3n
    J = genZeroN(numSheaths * numConfigParams);

    Vec tip_pos(3); // position of tip of robot
    Mat temp_matrix(3, 3);
    temp_matrix = genIdent3(); // placeholder matrix
    forwardKinematics(q, tip_pos, temp_matrix); // get global robot tip position
    
    for (int i = 0; i < numSheaths; ++i) {

        Vec prev_tip_pos(3);
        prev_tip_pos << long double(0.0), long double(0.0), long double(0.0); // tip position of previous segment chain
        Mat prev_rot = genIdent3(); // orientation matrix of tip of previous segment chain
        
        if (i > 0) {
            forwardKinematics(q.segment(0, i * numConfigParams), prev_tip_pos, prev_rot); // saves tip position and orientation of previous segment
        }

        Vec pos_ith_tip(3);
        pos_ith_tip << long double(0.0), long double(0.0), long double(0.0);                                                           // tip position of current segment
        Mat temp_matrix = genIdent3();                       // Not used in this function call

        forwardKinematics(q.segment(0, (i + 1) * (numConfigParams)), pos_ith_tip, temp_matrix); // saves tip position of current segment

        /** Jacobian formulation **/
        Mat M_i(6, 6);
        M_i << genIdent3(), long double(- 1.0) * skew(tip_pos - pos_ith_tip),
            genZero3(), genIdent3();

        Vec J_v_i_correspondingToTheta(3);
        Vec tempmultvec(3);
        tempmultvec << length(i) * (-1.0 * std::pow(theta(i), -2) + std::pow(theta(i), -2) * std::cos(theta(i)) + std::pow(theta(i), -1) * std::sin(theta(i))),
            0.0,
            length(i)* (-1.0 * std::pow(theta(i), -2) * std::sin(theta(i)) + std::pow(theta(i), -1) * std::cos(theta(i)));
        J_v_i_correspondingToTheta = prev_rot * getRotz_rad(phi(i)) * tempmultvec;

        Vec J_v_i_correspondingToPhi(3);
        tempmultvec << -1.0 * std::sin(phi(i)) * (1.0 - std::cos(theta(i))), std::cos(phi(i))* (1.0 - std::cos(theta(i))), 0.0;
        J_v_i_correspondingToPhi = prev_rot * (length(i) / theta(i)) * tempmultvec;

        Vec J_v_i_correspondingToLength(3);
        tempmultvec << std::pow(theta(i), -1) * (1.0 - std::cos(theta(i))), 0.0, std::pow(theta(i), -1)* std::sin(theta(i));
        J_v_i_correspondingToLength = prev_rot * getRotz_rad(phi(i)) * tempmultvec;

        Mat J_v_i(3, numConfigParams); // holds all upper components of the Jacobian
        J_v_i << J_v_i_correspondingToTheta, J_v_i_correspondingToPhi, J_v_i_correspondingToLength;

        Vec e2(3);
        e2 << long double(0.0), long double(1.0), long double(0.0);
        Vec e3(3);
        e3 << long double(0.0), long double(0.0), long double(1.0);
        Mat J_w_i(3, numConfigParams);
        J_w_i << prev_rot * getRotz_rad(phi(i)) * e2, prev_rot* e3, genZeros(3);

        Mat J_i(6, numConfigParams);
        J_i << J_v_i, J_w_i;

        // performs adjoint transformation M_i (6x6) on
        // each sub-Jacobian J_i (6x3) and places the
        // resulting 6x3 matrix into the J matrix
        J.block(0, i * numConfigParams, 6, numConfigParams) = M_i * J_i;
    }
    
    return J;
}

/*
    Calculates q_0_dot which is used for redundancy resolution by specifying
    the secondary task projected onto the nullspace of the Jacobian
*/
Vec generateNullSpaceObjectiveVector(Vec q,
    Vec jointMidPoints,
    Vec globalJointHighLimits,
    Vec globalJointLowLimits) {

    int numSheaths = q.rows() / numConfigParams;
    // joints in this case means the configuration parameters of the system
    Vec toReturn(numSheaths * numConfigParams); // holds secondary task that will be projected onto null space of Jacobian
    toReturn = genZeros(numSheaths * numConfigParams);

    Vec weights(numSheaths * numConfigParams);
    weights << 1.0, 1.0, 10.0, 1.0, 1.0, 10.0, 1.0, 1.0, 10.0;

    // calculate the secondary task for each configuration param
    // the secondary task for this system is to keep the joint (config param) values
    // as far away from the joint limits as possible -> therefore, keep the joint values
    // close to the midpoint of the joint limits
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (((i + 1) % 3) == 2) { // calculation for phi; not currently controlled
            toReturn(i) = long double(0.0);
        }
        else { // calculation for theta and length
            toReturn(i) = weights(i) * (-1.0 / (1.0 * (long double (numSheaths) * long double (numConfigParams)))) *
                (q(i) - jointMidPoints(i)) /
                ((globalJointHighLimits(i) - globalJointLowLimits(i)) *
                    (globalJointHighLimits(i) - globalJointLowLimits(i)));
        }
    }

    return toReturn;
}

/*
    Calculates and returns the joint velocities using the general inverse solution
    for differential kinematics
*/
Vec calculateJointVelocities(Mat J,
    Vec tipVelocityDesired, Vec nullspaceObjectiveVector) {
    double toler = 0.002; // used in rank estimation for pinv

    int numSheaths = J.cols() / numConfigParams;
    // calculate pseudoinverse of Jacobian
    Mat pinvJ = pinv(J, toler);
    // calculate image space
    Vec imageSpace = pinvJ * tipVelocityDesired;
    Mat ident = genIdentN(numSheaths * numConfigParams);
    // calculate nullspace
    Vec nullSpace = (ident - pinvJ * J) *
        redundancyResolutionOn *
        nullspaceObjectiveVector;
    //Vec nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
    //    redundancyResolutionOn *
    //    nullspaceObjectiveVector;
    // calculate joint velocities
    Vec qDot = imageSpace + nullSpace;

    bool skip = false;
    if ((nullSpace * nullSpace.transpose()).sum() == 0) {
        skip = false;
    }

    bool nullSpaceOverspeed = false;
    Vec nullspaceScaling(numSheaths * numConfigParams_);
    nullspaceScaling = genZeros(numSheaths * numConfigParams_); // scaling factor for downscaling the nullspace

    if (!skip) {
        for (int i = 0; i < (numSheaths * numConfigParams_); ++i) { // check every joint to see if it exceeds the designated max speed
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) { // calculate nullspace scaling factor
                long double nullSpaceElementIdealValue = sgn(qDot(i)) * maxJointSpeed(i) - imageSpace(i);
                nullspaceScaling(i) = nullSpace(i) / nullSpaceElementIdealValue;
                nullSpaceOverspeed = true;
            }
        }

        if (nullSpaceOverspeed) { // check that maximum scaling factor is greater than 0
            qDot = imageSpace + nullSpace / nullspaceScaling.maxCoeff(); // downscale nullspace
        }
    }
    else {
        qDot = imageSpace;
    }

    // Check again if new qDot is beyond the speed limit
    bool qDotOverspeed = false;
    Vec qDotScaling(numSheaths * numConfigParams);
    qDotScaling = genZeros(numSheaths * numConfigParams);
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000000000001)) {
            qDotOverspeed = true;
        }
    }

    if (qDotOverspeed) {
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000000000001)) {
                qDotScaling(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        qDot = qDot / qDotScaling.maxCoeff(); // apply recalculated scaling factors
    }

    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000001)) { // somehow still overspeeding
            std::cout << "overspeeding" << std::endl;
        }
    }

    return qDot;
}

Vec calculateJointVelocities(Mat J,
    Vec tipVelocityDesired, Vec nullspaceObjectiveVector, std::vector<double> lockedjoints, int exceededJointTracker) {
    double toler = 0.002; // used in rank estimation for pinv

    int numSheaths = J.cols() / numConfigParams;
    // calculate pseudoinverse of Jacobian
    Mat pinvJ = pinv(J, toler);
    // calculate image space
    Vec imageSpace = pinvJ * tipVelocityDesired;
    Mat ident = genIdentN(numSheaths * numConfigParams);
    // calculate nullspace
    Vec nullSpace = (ident - pinvJ * J) *
        redundancyResolutionOn *
        nullspaceObjectiveVector;
    //Vec nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
    //    redundancyResolutionOn *
    //    nullspaceObjectiveVector;
    // calculate joint velocities
    Vec qDot = imageSpace + nullSpace;

    bool skip = false;
    if ((nullSpace * nullSpace.transpose()).sum() == 0) {
        skip = true;
    }

    bool nullSpaceOverspeed = false;
    Vec nullspaceScaling(numSheaths * numConfigParams_);
    nullspaceScaling = genZeros(numSheaths * numConfigParams_); // scaling factor for downscaling the nullspace

    if (!skip) {
        for (int i = 0; i < (numSheaths * numConfigParams_); ++i) { // check every joint to see if it exceeds the designated max speed
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) { // calculate nullspace scaling factor
                long double nullSpaceElementIdealValue = sgn(qDot(i)) * maxJointSpeed(i) - imageSpace(i);
                nullspaceScaling(i) = nullSpace(i) / nullSpaceElementIdealValue;
                nullSpaceOverspeed = true;
            }
        }

        if (nullSpaceOverspeed) { // check that maximum scaling factor is greater than 0
            qDot = imageSpace + nullSpace / nullspaceScaling.maxCoeff(); // downscale nullspace
        }
    }
    else {
        qDot = imageSpace;
    }

    // Check again if new qDot is beyond the speed limit
    bool qDotOverspeed = false;
    Vec qDotScaling(numSheaths * numConfigParams);
    qDotScaling = genZeros(numSheaths * numConfigParams_);
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000000000001)) {
            qDotOverspeed = true;
        }
    }

    if (qDotOverspeed) {
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000000000001)) {
                qDotScaling(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        qDot = qDot / qDotScaling.maxCoeff(); // apply recalculated scaling factors
    }

    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + long double(0.0000000001)) { // somehow still overspeeding
            std::cout << "overspeeding" << std::endl;
        }
    }

    return qDot;
}


Vec inverseKinematics(Vec &q, Vec &qDot,
    Vec &displacementCalc,
    Vec tipVelocityDesired, int mode, int& outOfWorkspace) {
    // these must be hardcoded
    Vec cushion(3); // deadband around joint limits
    cushion << 0.01, 0.0, 0.5; // order is theta, phi, length with units of rad, rad, mm

    int numSheaths = q.rows() / numConfigParams;
    // these must be hardcoded
    Vec globalJointLowLimits(numConfigParams * numSheaths); // minimum allowable lengths/angles of the joints
    globalJointLowLimits << long double(0.0000000000000001), long double(-1000.0 * M_PI), lengthLowLimit, -q(5) / minAllowedRadii, long double(-1000.0 * M_PI), lengthLowLimit, -q(8) / minAllowedRadii, long double(-1000.0 * M_PI), lengthLowLimit;

    // these must be hardcoded
    Vec globalJointHighLimits(numConfigParams * numSheaths); // maximum allowable lengths/angles of the joints
    globalJointHighLimits << q(2) / minAllowedRadii, long double(1000.0 * M_PI), lengthHighLimit, q(5) / minAllowedRadii, long double(1000.0 * M_PI), lengthHighLimit, q(8) / minAllowedRadii, long double(1000.0 * M_PI), lengthHighLimit;
    
    // perform differential kinematics calculation to obtain Jacobian
    Mat J_normal = generateJacobian(q);

    Mat J_use; // Jacobian that will be used for robot control
    Vec tipVelocityDesired_use; // desired input that will be used for robot control
    std::vector<double> locked_joints; // non-actuated joints

    switch (mode) {
    case 0: // only position (3 DOF) with 5 DOF (no phi1 and sheath3)
        {
            locked_joints.push_back(1);
            locked_joints.push_back(6);
            locked_joints.push_back(7);
            locked_joints.push_back(8);
            // only grabs the top half of the Jacobian ie linear components
            J_use = J_normal.block(0, 0, 3, numConfigParams_ * numSheaths);
            // only grabs the top half of desired input ie x, y, z
            tipVelocityDesired_use = tipVelocityDesired.head(3);
            break;
        }
    case 1: // position and orientation with no roll (5 DOF) with 6 DOF robot (no phi1, theta3, and phi3)
        {
            locked_joints.push_back(1);
            locked_joints.push_back(6);
            locked_joints.push_back(7);
            // J_use = J_normal.block(0, 0, 5, numConfigParams_ * numSheaths);
            // tipVelocityDesired_use = tipVelocityDesired.head(5);
            Mat temp1(4, numConfigParams_ * numSheaths);
            temp1 << J_normal.row(0), J_normal.row(1), J_normal.row(2), J_normal.row(4);
            J_use = temp1;

            Vec temp2(4);
            temp2 << tipVelocityDesired(0), tipVelocityDesired(1), tipVelocityDesired(2), tipVelocityDesired(4);
            tipVelocityDesired_use = temp2;
            break;
        }
    default: // position and orientation (6 DOF) with 7 DOF robot (no phi1 and theta3)
        {
            locked_joints.push_back(1);
            locked_joints.push_back(6);
            J_use = J_normal;
            tipVelocityDesired_use = tipVelocityDesired;
        }
    }

    // calculate secondary task for null space manipulation (redundancy resolution)
    Vec nullspaceObjectiveVector = generateNullSpaceObjectiveVector(q,
        jointMidPoints,
        globalJointHighLimits,
        globalJointLowLimits);

    for (int i : locked_joints) {
        J_use.col(i) = genZeros(J_use.rows());
        nullspaceObjectiveVector(i) = 0.0;

        switch (i % numConfigParams) {
        case 1: // phi, units in rad
            if (i == 1) {
                globalJointLowLimits(i) = long double(M_PI);
                globalJointHighLimits(i) = long double(M_PI);
            }
            else {
                globalJointLowLimits(i) = 0.0;
                globalJointHighLimits(i) = 0.0;
            }
            break;
        case 2: // length, units in mm
            globalJointLowLimits(i) = 1.0;
            globalJointHighLimits(i) = 1.0;
            break;
        default: // theta, units in rad
            globalJointLowLimits(i) = 0.01;
            globalJointHighLimits(i) = 0.01;
        }
    }

    // calculate joint velocities
    qDot = calculateJointVelocities(J_use,
        tipVelocityDesired_use,
        nullspaceObjectiveVector,
        locked_joints,
        -1);

    Vec qNew = q + qDot * deltaT; // recalculate config params

    int loop_counter = 0;
    Vec joints_exceeded(numSheaths * numConfigParams);
    joints_exceeded = genZeros(numSheaths * numConfigParams_);

    if (applyJointLimits) {
        bool finishedChecking = false;
        int numJoints = numSheaths * numConfigParams;
        Vec limitHighJoint = globalJointHighLimits;
        Vec limitLowJoint = globalJointLowLimits;

        while (!finishedChecking) {
            // calculate joint limits dynamically
            for (int i = 0; i < numSheaths; ++i) {
                if (i > 0) {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numConfigParams), -1.0 * (q(i * numConfigParams + 2) / minAllowedRadii)); // theta
                }
                else {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numConfigParams), long double(0.0)); // theta
                }
                limitLowJoint(i * numConfigParams + 2) = std::max(globalJointLowLimits(i * numConfigParams + 2), (q(i * numConfigParams) * minAllowedRadii)); // length

                limitHighJoint(i * numConfigParams) = std::min(globalJointHighLimits(i * numConfigParams), (q(i * numConfigParams + 2) / minAllowedRadii)); // theta
                limitHighJoint(i * numConfigParams + 2) = globalJointHighLimits(i * numConfigParams + 2); // length
            }

            int exceededJoint = 10;
            // 
            // find and save first joint that exceeds previously calculated joint limits
            for (int j = 0; j < numJoints; ++j) {
                if ((j % numConfigParams) != 1 &&
                    (find_element(locked_joints, j) == -1) &&
                    joints_exceeded(j) == 0) {
                    if (qNew(j) < (limitLowJoint(j) + cushion((j % numConfigParams))) ||
                        qNew(j) > (limitHighJoint(j) - cushion(j % numConfigParams))) {
                        exceededJoint = j;
                        if (qNew(j) > limitLowJoint(j) || qNew(j) < limitHighJoint(j)) {
                            joints_exceeded(j) = 1;
                        }
                        break;
                    }
                }
            }

            if (exceededJoint != 10) { // if a joint has exceeded the limits
                //unemployed++;
                J_use.col(exceededJoint) = genZeros(J_use.rows()); // wipe the column in the Jacobian of the current param
                nullspaceObjectiveVector(exceededJoint) = 0.0; // wipe q_0_dot

                // calculate joint velocities
                qDot = calculateJointVelocities(J_use,
                    tipVelocityDesired_use,
                    nullspaceObjectiveVector,
                    locked_joints,
                    exceededJoint);

                qNew = q + qDot * deltaT; // recalculate config params
            }
            else {
                finishedChecking = true;
            }
            loop_counter++;
        }
    }

    // update tendon displacements
    //displacementCalc = jointSpaceToActuationSpace(qNew);

    outOfWorkspace = (loop_counter > 1);
    return qNew;
}

std::pair<Vec, Vec> calculate_circleTrajectory_position_velocity_f(Vec start_point,
    long double elapsed_time,
    long double tip_speed,
    int out_of_workspace,
    Vec pos_cur) {
    long double radius = 20.0;
    Vec euler_angles(3);
    euler_angles << 0.0, M_PI / 4, 0.0;
    Vec circle_center_point(3);
    circle_center_point << start_point(0) - radius,
        start_point(1),
        start_point(2);

    //std::cout << "euler_angles: " << euler_angles << std::endl;
    //std::cout << "circle_center_point: " << circle_center_point << std::endl;

    // Calculate position 
    long double omega = tip_speed / radius;
    long double theta = omega * elapsed_time;

    //std::cout << "omega: " << omega << std::endl;
    //std::cout << "theta: " << theta << std::endl;

    long double x = radius * std::cos(theta) + circle_center_point(0);
    long double y = radius * std::sin(theta) + circle_center_point(1);
    long double z = circle_center_point(2);

    //std::cout << "x: " << x << std::endl;
    //std::cout << "y: " << y << std::endl;
    //std::cout << "z: " << z << std::endl;

    // Calculate velocity
    long double vx = -tip_speed * std::sin(theta);
    long double vy = tip_speed * std::cos(theta);
    long double vz = 0.0;

    // Rotation Matrices
    Mat Rz = getRotz_rad(euler_angles(2));

    Mat Ry = getRoty_rad(euler_angles(1));

    Mat Rx = getRotx_rad(euler_angles(0));

    //std::cout << "Rx: " << Rx << std::endl;
    //std::cout << "Ry: " << Ry << std::endl;
    //std::cout << "Rz: " << Rz << std::endl;

    // Final Rotation Matrix
    Mat rotation_matrix = Rz * Ry * Rx;

    // Rotate Velocity
    Vec temp(3);
    temp << vx, vy, vz;
    Vec velocity = rotation_matrix * temp;

    // Output position
    temp << x, y, z;
    Vec position = rotation_matrix * (temp - start_point) + start_point;

    if (out_of_workspace == 1) {
        velocity(0) = (position(0) - pos_cur(0)) / deltaT;
        velocity(1) = (position(1) - pos_cur(1)) / deltaT;
        velocity(2) = (position(2) - pos_cur(2)) / deltaT;
    }

    // Return results using std::pair
    return std::make_pair(position, velocity);
}

double solve_equation(long double f, long double a, long double b, long double c) {
    // Initial guess (you may need to adjust this based on your specific problem)
    long double th = -10.0;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        long double func = a + b * std::atan(th) + c * th - f;
        long double derivative = b / (1 + th * th) + c;

        long double delta = func / derivative;
        th -= delta;

        if (std::abs(delta) < LDBL_EPSILON) {
            return th;
        }
    }

    // If we reach here, the method didn't converge
    std::cerr << "Solution did not converge" << std::endl;
    return NAN;
}

double f1(double th, double phi, double a1, double b1, double c1) {
    return (a1 + b1 * std::atan(th) + c1) * std::cos(phi);
}

double f2(double th, double phi, double a2, double b2, double c2) {
    return (a2 + b2 * std::atan(th) + c2) * std::sin(phi);
}

std::vector<double> solveEquations(double f1_target, double f2_target,
    double a1, double b1, double c1,
    double a2, double b2, double c2) {
    double th = 26.0, phi = deg2rad(56.0); // Initial guesses
    double epsilon = 1e-6; // Tolerance for convergence
    int maxIterations = 1000;

    for (int i = 0; i < maxIterations; ++i) {
        // Calculate function values
        double F1 = f1(th, phi, a1, b1, c1) - f1_target;
        double F2 = f2(th, phi, a2, b2, c2) - f2_target;

        // Calculate Jacobian elements
        double J11 = b1 * std::cos(phi) / (1 + th * th);
        double J12 = -(a1 + b1 * std::atan(th) + c1) * std::sin(phi);
        double J21 = b2 * std::sin(phi) / (1 + th * th);
        double J22 = (a2 + b2 * std::atan(th) + c2) * std::cos(phi);

        // Calculate determinant of Jacobian
        double det = J11 * J22 - J12 * J21;

        // Calculate updates
        double delta_th = (-F1 * J22 + F2 * J12) / det;
        double delta_phi = (-F2 * J11 + F1 * J21) / det;

        // Update th and phi
        th += delta_th;
        phi += delta_phi;

        // Check for convergence
        if (std::abs(delta_th) < epsilon && std::abs(delta_phi) < epsilon) {
            return { th, phi };
        }
    }

    // If no convergence, return empty vector
    return {};
}

Eigen::Vector3d rotationMatrixToEulerAngles(Eigen::Matrix3d R, double epsilon = 1e-6) {
    // Check if the rotation matrix is valid

    double x, y, z;

    // Check for gimbal lock
    if (std::abs(R(0,2)) > 1.0 - epsilon) {
        // Gimbal lock case
        z = 0.0;
        if (R(0,2) < 0) {
            y = PI / 2.0;
            x = std::atan2(R(1,0), R(2,0));
        }
        else {
            y = -PI / 2.0;
            x = std::atan2(-R(1,0), -R(2,0));
        }
    }
    else {
        // General case
        y = std::asin(-R(0,2));
        double cos_y = std::cos(y);

        x = std::atan2(R(1,2) / cos_y, R(2,2) / cos_y);
        z = std::atan2(R(0,1) / cos_y, R(0,0) / cos_y);
    }

    // Normalize angles to be between -PI and PI
    x = std::fmod(x + M_PI, 2 * M_PI) - M_PI;
    y = std::fmod(y + M_PI, 2 * M_PI) - M_PI;
    z = std::fmod(z + M_PI, 2 * M_PI) - M_PI;

    return { y, x, z }; // Return in YXZ order
}


Eigen::Vector3d customRotToAng(Eigen::Matrix3d rot, double ry_p, double rx_p, double rz_p) {
    int i = 1;
    int j = 0;
    int k = 2;

    Eigen::Vector3d toReturn = Eigen::Vector3d::Zero();

    toReturn(0) = std::atan2(rot(j, k), rot(k, k));
    auto c2 = Eigen::Vector2d(rot(i, i), rot(i, j)).norm();

    if (toReturn(0) < 0.0) {
        //toReturn(0) += M_PI;
        toReturn(1) = std::atan2(-1.0 * rot(i, k), -1.0 * c2);
    } else {
        toReturn(1) = std::atan2(-1.0 * rot(i, k), c2);
    }

    auto s1 = std::sin(toReturn(0));
    auto c1 = std::cos(toReturn(0));
    toReturn(2) = std::atan2(s1 * rot(k, i) - c1 * rot(j, i), c1 * rot(j, j) - s1 * rot(k, j));

    return toReturn;
}

Eigen::Vector3d rotationToEuler_yxz(Eigen::Matrix3d toConvert, double prev_rx, double prev_ry, double prev_rz) {
    Eigen::Vector3d toReturn;
    double x = std::asin(-1.0 * toConvert(1, 2));
    double y = std::atan2(toConvert(0, 2) / std::cos(x), toConvert(2, 2) / std::cos(x));
    double z = std::atan2(toConvert(1, 0) / std::cos(x), toConvert(1, 1) / std::cos(x));

    if (std::abs(y - prev_ry) > 0.2) {
        x = (x > 0) ? (M_PI - x) : (-1.0 * x - M_PI);

        y = std::atan2(toConvert(0, 2) / std::cos(x), toConvert(2, 2) / std::cos(x));
        z = std::atan2(toConvert(1, 0) / std::cos(x), toConvert(1, 1) / std::cos(x));
    }

    toReturn << y, x, z;
    return toReturn;
}

int main()
{
    fesetround(FE_TONEAREST); // Round to nearest (even) like MATLAB

    theta_init_ << deg2rad(long double(25.0)), deg2rad(long double(25.0)), long double(0.00005); // initial bending angles of each sheath
    phi_init_ << long double( 180.0 * M_PI / 180.0), long double( - 90.0 * M_PI / 180.0), long double(0.0 * M_PI / 180.0); // initial bending plane offsets of each sheath
    length_init_ << long double(10.0), long double(10.0), long double(10.0); // initial length of the steerable section of each sheath in mm
    delta_ << long double(2.3), long double(1.8), long double(1.0); // radius (in meters) of cross sectional circle centered at sheath centerline and passes through the tendon locations

    Vec q(numConfigParams * numTubes);
    Vec qDot(numConfigParams * numTubes);

    q_sim << long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0);
    qDot_sim << long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0);
    displacement_sim << long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0), long double(0.0);

    jointMidPoints << long double(M_PI / 4.0), long double(0.0), long double(25.0), long double(0.0), long double(0.0), long double(25.0), long double(0.0), long double(0.0), long double(25.0);

    maxJointSpeed << long double(M_PI / 6.0), long double(1000.0 * M_PI), long double(3.0),
        long double(M_PI / 6.0), long double(1000.0 * M_PI), long double(3.0),
        long double(M_PI / 6.0), long double(1000.0 * M_PI), long double(3.0);

    Mat tendonAngOffsets(numTotalTubes, numTendonsPerTube);
    tendonAngOffsets.row(0) << long double(-1.0 * M_PI), long double(-1.0 * 3.0 * M_PI / 2.0); // outer sheath
    tendonAngOffsets.row(1) << long double(-1.0 * M_PI / 4.0), long double(-1.0 * 3.0 * M_PI / 4.0); // intermediate sheath
    tendonAngOffsets.row(2) << long double(-1.0 * 0.0), long double(-1.0 * M_PI / 2.0); // TODO: inner sheath (currently unused)

    for (int i = 0; i < numTotalTubes; ++i) {
        tendonGeometry[i] << std::cos(tendonAngOffsets(i, 0)),
            std::cos(tendonAngOffsets(i, 1)),
            std::sin(tendonAngOffsets(i, 0)),
            std::sin(tendonAngOffsets(i, 1));

        invTendonGeometry[i] = tendonGeometry[i].inverse();
    }

    for (int i = 0; i < numTubes; ++i) {
        q_sim.segment((i * numConfigParams), numConfigParams) << theta_init_(i), phi_init_(i), length_init_(i);
    }

    std::vector<MedianFilter> medfilt;
    for (int i = 0; i < 6; ++i) {
        medfilt.push_back(MedianFilter(100));
    }

    std::cout << q_sim << std::endl;
    Mat dummyRotMM_i_wrtG1(3, 3);
    dummyRotMM_i_wrtG1 = genIdent3();

    //Eigen::VectorXd q_test(9);
    //q_test << deg2rad(67), deg2rad(54), 20, deg2rad(67), deg2rad(54), 20, deg2rad(67), deg2rad(54), 20;

    forwardKinematics(q_sim, pTip_init_sim, dummyRotMM_i_wrtG1);
    //std::cout << "ptip: " << pTip_init_sim << std::endl;
    //std::cout << "rot: " << dummyRotMM_i_wrtG1 << std::endl;

    //std::cout << generateJacobian(q_sim) << std::endl;

    //std::vector<Eigen::VectorXd> printqsims;
    //std::vector<Eigen::VectorXd> printqdotsims;

    automated_velocity << long double(0.0), long double(0.0), long double(0.0);

    for (int j = 0; j < 10000; ++j) {
        auto time = j * deltaT;
        Mat dummyRotMM_i_wrtG(3, 3);
        dummyRotMM_i_wrtG = genIdent3();
        
        forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);

        std::pair<Vec, Vec> result =
            calculate_circleTrajectory_position_velocity_f(pTip_init_sim, time, long double(0.25), outOfWorkspace, pTip_sim);
        automated_velocity = result.second;

        Vec tipVelocityDesired(6); // command from manual control
        tipVelocityDesired << automated_velocity(0), automated_velocity(1), automated_velocity(2), long double(0.0), long double(0.1745 * (2.0 * M_PI / 500.0) * std::cos(j * ((2.0 * M_PI) / 500.0))), long double (0.0);

        Vec displacementCalc(nOfDrivenMotors);
        displacementCalc = genZeros(nOfDrivenMotors); // motor displacement

        //if (j >= 15000 && j <= 15500) {
        //    std::cout << j << std::endl;
        //    std::cout << "qsim: " << q_sim << std::endl;
        //    std::cout << "qDot: " << qDot_sim << std::endl;
        //    std::cout << "***********************" << std::endl;
        //}

        //writeVectorToCSV("C:/Users/ch256744/Downloads/q_output_precise.csv", q_sim);
        //writeVectorToCSV("C:/Users/ch256744/Downloads/qdot_output_precise.csv", qDot_sim);
        //writeVectorToCSV("C:/Users/ch256744/Downloads/tipspeed_output_precise.csv", tipVelocityDesired);
        

        q_sim = inverseKinematics(q_sim, qDot_sim, displacement_sim, tipVelocityDesired, operatingMode, outOfWorkspace);
        Mat jacobian = generateJacobian(q_sim);


        forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);
        writeVectorToCSV("C:/Users/ch256744/BCH Dropbox/Phillip Tran/ExtendedICRAPaper/AbdulsCode/MatlabCode/endtip_longdouble1021_3.csv", pTip_sim);
        unemployed += outOfWorkspace;
    }

    //q_sim << 0.657331119677952, 3.14159265358979, 46.1281520431707, 0.742844152344368, -1.74947466638556, 49.4991203380709, 5.00000000000000e-05, 0.0, 30.0;
    //Eigen::VectorXd tipVelocityDesired(6);
    //tipVelocityDesired << -0.225868151611763, -0.267126373899629, 0.202975859771470, 0.0, 0.0, 0.0;

    //q_sim = inverseKinematics(q_sim, qDot_sim, displacement_sim, tipVelocityDesired, operatingMode, outOfWorkspace);

    std::cout << "new ptip: " << pTip_sim << std::endl;
    std::cout << "qsim: " << q_sim << std::endl;
    std::cout << "unemployed: " << unemployed << std::endl;
    std::cout << "qDot: " << qDot_sim << std::endl;

    std::cout << "test: " << (2 > 1) << std::endl;

    //Eigen::MatrixXd j(3, 9);
    //j << -43.376, 0, 0, 22.851, -25.555, 0, 0, 0, 0,
    //    -6.3036e-13, 0, 0, 14.115, 1.35, 0, 0, 0, 0,
    //    -39.336, 0, 0, -18.3, -30.869, 0, 0, 0, 0;
    //std::cout << "j " << j << std::endl;
    //std::cout << "jinv" << pinv(j, 0.002) << std::endl;

    Eigen::MatrixXd pinv(9,3);
    pinv << 0.483791, -2.47528, -0.561424, 0, 0, 0, 0, 0, 0, 0.186841, -0.821548, -0.205854, -0.597136, 2.98518, 0.657902, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd j(3, 9);
    j << -38.7923, 0, 0, 24.8744, -25.3205, 0, 0, 0, 0, 4.75069e-15, 0, 0, 8.88868, 2.78123, 0, 0, 0, 0, -35.2094, 0, 0, -17.7548, -34.0815, 0, 0, 0, 0;
    std::cout << "h: " << pinv * j << std::endl;
    Eigen::VectorXd ns(9);
    ns << -0.000956319, 0, 0, -0.00273355, 0, 0, 0, 0, 0;
    Eigen::MatrixXd p = Eigen::MatrixXd::Identity(9,9);
    std::cout << "p: " << p - (pinv*j) << std::endl;
    std::cout << "g: " << (p - (pinv * j)) * ns << std::endl;

    std::cout << "vec size: " << ns.size() << std::endl;
};

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
