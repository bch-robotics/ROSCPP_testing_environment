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
#include <random>

typedef Eigen::Matrix< long double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > Mat;

class SlidingModeFilter {
private:
    double lambda;
    double k1;
    double k2;
    double T;
    double x1;
    double x2;

    double sign(double value) {
        return (value > 0) ? 1.0 : ((value < 0) ? -1.0 : 0.0);
    }

public:
    SlidingModeFilter(double lambda, double k1, double k2, double T)
        : lambda(lambda), k1(k1), k2(k2), T(T), x1(0.0), x2(0.0) {}

    void reset() {
        x1 = 0.0;
        x2 = 0.0;
    }

    double filter(double input) {
        double e = input - x1;
        double s = e + lambda * x2;

        x1 = x1 + T * x2;
        x2 = x2 - T * (k1 * sign(s) + k2 * s);

        return x1;
    }

    std::vector<double> filterSignal(const std::vector<double>& inputSignal) {
        std::vector<double> outputSignal;
        outputSignal.reserve(inputSignal.size());

        for (const auto& input : inputSignal) {
            outputSignal.push_back(filter(input));
        }

        return outputSignal;
    }
};

class MedianFilter {
private:
    std::vector<double> buffer;
    size_t windowSize;

public:
    MedianFilter(size_t size) : windowSize(size) {
        if (size % 2 == 0 || size < 3) {
            throw std::invalid_argument("Window size must be an odd number greater than or equal to 3");
        }
        buffer.reserve(size);
    }

    double filter(double input) {
        if (buffer.size() < windowSize) {
            buffer.push_back(input);
        }
        else {
            buffer.erase(buffer.begin());
            buffer.push_back(input);
        }

        std::vector<double> sortedBuffer = buffer;
        std::sort(sortedBuffer.begin(), sortedBuffer.end());

        size_t medianIndex = sortedBuffer.size() / 2;
        return sortedBuffer[medianIndex];
    }

    void reset() {
        buffer.clear();
    }

    size_t getWindowSize() const {
        return windowSize;
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
const double EPSILON = 1e-6;
const int MAX_ITERATIONS = 100;

Eigen::Vector3d automated_velocity = Eigen::Vector3d::Zero();
Eigen::Matrix<double, 2, 2, Eigen::RowMajor> tendonGeometry[3]; // creates a 3x2x2 matrix ie 3 2x2 matrices
Eigen::Matrix<double, 2, 2, Eigen::RowMajor> invTendonGeometry[3]; // creates a 3x2x2 matrix ie 3 2x2 matrices

Eigen::VectorXd theta_init_(numConfigParams); // initial bending angles of each sheath
Eigen::VectorXd phi_init_(numConfigParams); // initial bending plane offsets of each sheath
Eigen::VectorXd length_init_(numConfigParams); // initial length of the steerable section of each sheath in mm
Eigen::VectorXd delta_(numConfigParams); // radius (in meters) of cross sectional circle centered at sheath centerline and passes through the tendon locations
double lengthLowLimit = 1.0;  // in mm
double lengthHighLimit = 50.0; // in mm
double minAllowedRadii = 12.0; // in mm, calculated using length = radius * theta (in radians) where theta is 90 degrees
Eigen::VectorXd maxJointSpeed = Eigen::VectorXd::Zero(numTubes * numConfigParams);
Eigen::VectorXd jointMidPoints = Eigen::VectorXd::Zero(numTubes * numConfigParams);
int outOfWorkspace = 0;
double previous_time_ = 0.0;
double deltaT = 0.1;
double redundancyResolutionOn = 1.0;
bool useOnlyPositionJacobian = true;
bool applyJointLimits = true;
/*  operating modes
    0: only position (3 DOF) with 5 DOF (no phi1 and sheath3),
    1: position and orientation with no roll (5 DOF) with 6 DOF robot (no phi1, theta3, and phi3)
    2: position and orientation (6 DOF) with 7 DOF robot (no phi1 and theta3).
*/
int operatingMode = 1;

bool exp_start_time_saved_ = false;
double exp_start_time_ = 0.0;
int unemployed = 0;

// FOR DEBUGGING AND SIM 
Eigen::VectorXd q_sim = Eigen::VectorXd::Zero(numTubes * numConfigParams);
Eigen::VectorXd qDot_sim = Eigen::VectorXd::Zero(numTubes * numConfigParams);
Eigen::Vector3d pTip_init_sim;
Eigen::Vector3d pTip_sim;
Eigen::VectorXd displacement_sim = Eigen::VectorXd::Zero(nOfDrivenMotors);

template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Define the sign function
double sign(double value) {
    if (value > 0) return 1.0;
    if (value < 0) return -1.0;
    return 0.0;
}

// Implement the 3rd order sliding mode filter
std::vector<double> slidingModeFilter(const std::vector<double>& input, double lambda, double gain, double dt) {
    size_t n = input.size();
    std::vector<double> output(n, 0.0);

    // Initial conditions
    double x1 = 0.0, x2 = 0.0, x3 = 0.0;

    for (size_t i = 1; i < n; ++i) {
        // Calculate the sliding variable
        double sigma = input[i] - x1;

        // Update the states using a simple Euler method
        x1 += dt * x2;
        x2 += dt * x3;
        x3 += dt * (-gain * sign(sigma + lambda * x2 + lambda * lambda * x3));

        // Store the output
        output[i] = x1;
    }

    return output;
}

Eigen::Matrix3d skew(const Eigen::Vector3d& w)
{
    if (3 != w.size())
    {
        throw std::invalid_argument("Vector must be 3x1");
    }

    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();

    R(0, 1) = -1.0 * w(2);
    R(0, 2) = w(1);
    R(1, 2) = -1.0 * w(0);

    R(1, 0) = w(2);
    R(2, 0) = -1.0 * w(1);
    R(2, 1) = w(0);

    return R;
}

Eigen::MatrixXd pinv(Eigen::MatrixXd A, double tol)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = svd.singularValues();

    ////std::cout << "single ladies" << s << std::endl;

    if (!s.isZero(0))
    {
        // int r1 = (s.array() > tol).count() + 1; // since the last s is too small and this +1 allows it to exist, its inverse down the line creates big values. This is to ensure that r1 is at least 1 but it is now done as follows (but matlab's pinv if you look at it carefully does not do it. although there is +1):
        int r1 = (s.array() > tol).count(); // estimate effective rank
        r1 = std::max(r1, 1);               // Ensure that r1 is at least 1 // matlab does not ensure this

        Eigen::MatrixXd U = svd.matrixU();
        Eigen::MatrixXd V = svd.matrixV();

        // V.rightCols(V.cols() - r1).setZero();
        V.conservativeResize(V.rows(), r1);
        // U.rightCols(U.cols() - r1).setZero();
        U.conservativeResize(U.rows(), r1);
        s.conservativeResize(r1);

        s = s.array().inverse();
        return V * s.asDiagonal() * U.transpose();
    }
    else
    {
        return Eigen::MatrixXd::Constant(A.cols(), A.rows(), std::numeric_limits<double>::quiet_NaN());
    }

    //auto toReturn = A.completeOrthogonalDecomposition();
    //toReturn.setThreshold(tol); // sets the threshold for which values should be considered to be 0
    //return Eigen::MatrixXd(toReturn.pseudoInverse());
}

double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

double rad2deg(double rad) {
    return rad * 180.0 / M_PI;
}

Eigen::Vector3d rad2deg(Eigen::Vector3d rad) {
    Eigen::Vector3d toReturn;
    toReturn << rad2deg(rad(0)), rad2deg(rad(1)), rad2deg(rad(2));
    return toReturn;
}

Eigen::Matrix3d getRotx(double ang) {
    Eigen::Matrix3d toReturn;
    toReturn << 1.0, 0.0, 0.0,
        0.0, std::cos(deg2rad(ang)), -1.0 * std::sin(deg2rad(ang)),
        0.0, std::sin(deg2rad(ang)), std::cos(deg2rad(ang));

    return toReturn;
}

Eigen::Matrix3d getRoty(double ang) {
    Eigen::Matrix3d toReturn;
    toReturn << std::cos(deg2rad(ang)), 0.0, std::sin(deg2rad(ang)),
        0.0, 1.0, 0.0,
        -1.0 * std::sin(deg2rad(ang)), 0.0, std::cos(deg2rad(ang));

    return toReturn;
}

Eigen::Matrix3d getRotz(double ang) {
    Eigen::Matrix3d toReturn;
    toReturn << std::cos(deg2rad(ang)), -1.0 * std::sin(deg2rad(ang)), 0.0,
        std::sin(deg2rad(ang)), std::cos(deg2rad(ang)), 0.0,
        0.0, 0.0, 1.0;

    return toReturn;
}

Eigen::Matrix3d getRotx_rad(double ang_rad) {
    return getRotx(rad2deg(ang_rad));
}

Eigen::Matrix3d getRoty_rad(double ang_rad) {
    return getRoty(rad2deg(ang_rad));
}

Eigen::Matrix3d getRotz_rad(double ang_rad) {
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

void readCSVDub(const std::string& filename, std::vector< std::vector<double> >& toReturn) {
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
                //std::cout << value << std::endl;
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
            toReturn.push_back(row);
        }
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

void writeCSV(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    for (const auto& value : data) {
        file << std::fixed << std::setprecision(6) << value << "\n";
    }

    file.close();
}

void writeVectorToCSV(const std::string& filename, const Eigen::VectorXd& vector) {
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
void forwardKinematics(Eigen::VectorXd q,
    Eigen::Vector3d& end_tip_pos,
    Eigen::Matrix3d& rotationMatrix) {

    int numSheaths = q.rows() / numConfigParams;

    Eigen::VectorXd theta(numSheaths);
    Eigen::VectorXd phi(numSheaths);
    Eigen::VectorXd length(numSheaths);

    // store the config space params in separate vectors
    for (int i = 0; i < numSheaths; ++i) {
        theta(i) = (q(i * numConfigParams) == 0.0) ? std::numeric_limits<double>::epsilon() : q(i * numConfigParams);
        phi(i) = q(i * numConfigParams + 1);
        length(i) = q(i * numConfigParams + 2);
    }

    // initialize rotation matrix and end tip position vector
    rotationMatrix = Eigen::Matrix3d::Identity();
    end_tip_pos.setZero();

    for (int k = 0; k < numSheaths; ++k) {
        rotationMatrix = Eigen::Matrix3d::Identity();

        for (int j = 0; j < (k + 1); ++j) {

            if (j < k) { // calculate rotation matrix if not at final sheath
                rotationMatrix = rotationMatrix * (getRotz_rad(phi(j))) * (getRoty_rad(theta(j)));
            }
            else { // calculate end tip position of current sheath wrt previous frame
                Eigen::Vector3d tip_pos_rel;
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
Eigen::MatrixXd generateJacobian(Eigen::VectorXd q) {

    int numSheaths = q.rows() / numConfigParams;

    Eigen::VectorXd theta(numSheaths);
    Eigen::VectorXd phi(numSheaths);
    Eigen::VectorXd length(numSheaths);

    // separate and store config params
    for (int i = 0; i < numSheaths; ++i) {
        theta(i) = q(i * numConfigParams);
        phi(i) = q(i * numConfigParams + 1);
        length(i) = q(i * numConfigParams + 2);
    }

    Eigen::MatrixXd J(6, numSheaths * numConfigParams); // Jacobian matrix, 6x3n
    J.setZero();

    Eigen::Vector3d tip_pos; // position of tip of robot
    Eigen::Matrix3d temp_matrix = Eigen::Matrix3d::Identity(); // placeholder matrix
    forwardKinematics(q, tip_pos, temp_matrix); // get global robot tip position
    
    for (int i = 0; i < numSheaths; ++i) {

        Eigen::Vector3d prev_tip_pos = { 0.0, 0.0, 0.0 }; // tip position of previous segment chain
        Eigen::Matrix3d prev_rot = Eigen::Matrix3d::Identity(); // orientation matrix of tip of previous segment chain
        
        if (i > 0) {
            forwardKinematics(q.segment(0, i * numConfigParams), prev_tip_pos, prev_rot); // saves tip position and orientation of previous segment
        }

        Eigen::Vector3d pos_ith_tip = { 0.0, 0.0, 0.0 };                                                           // tip position of current segment
        Eigen::Matrix3d temp_matrix = Eigen::Matrix3d::Identity();                       // Not used in this function call

        forwardKinematics(q.segment(0, (i + 1) * (numConfigParams)), pos_ith_tip, temp_matrix); // saves tip position of current segment

        /** Jacobian formulation **/
        Eigen::MatrixXd M_i(6, 6);
        M_i << Eigen::Matrix3d::Identity(), -1 * skew(tip_pos - pos_ith_tip),
            Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

        Eigen::Vector3d J_v_i_correspondingToTheta = prev_rot * getRotz_rad(phi(i)) *
            Eigen::Vector3d(
                length(i) * (-1.0 * std::pow(theta(i), -2) + std::pow(theta(i), -2) * std::cos(theta(i)) + std::pow(theta(i), -1) * std::sin(theta(i))),
                0.0,
                length(i) * (-1.0 * std::pow(theta(i), -2) * std::sin(theta(i)) + std::pow(theta(i), -1) * std::cos(theta(i))));

        Eigen::Vector3d J_v_i_correspondingToPhi = prev_rot * (length(i) / theta(i)) * Eigen::Vector3d(-1.0 * std::sin(phi(i)) * (1.0 - std::cos(theta(i))), std::cos(phi(i)) * (1.0 - std::cos(theta(i))), 0.0);

        Eigen::Vector3d J_v_i_correspondingToLength = prev_rot * getRotz_rad(phi(i)) * Eigen::Vector3d(std::pow(theta(i), -1) * (1.0 - std::cos(theta(i))), 0.0, std::pow(theta(i), -1) * std::sin(theta(i)));

        Eigen::MatrixXd J_v_i(3, numConfigParams); // holds all upper components of the Jacobian
        J_v_i << J_v_i_correspondingToTheta, J_v_i_correspondingToPhi, J_v_i_correspondingToLength;

        Eigen::Vector3d e2(0.0, 1.0, 0.0);
        Eigen::Vector3d e3(0.0, 0.0, 1.0);
        Eigen::MatrixXd J_w_i(3, numConfigParams);
        J_w_i << prev_rot * getRotz_rad(phi(i)) * e2, prev_rot* e3, Eigen::Vector3d::Zero();

        Eigen::MatrixXd J_i(6, numConfigParams);
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
Eigen::VectorXd generateNullSpaceObjectiveVector(Eigen::VectorXd q,
    Eigen::VectorXd jointMidPoints,
    Eigen::VectorXd globalJointHighLimits,
    Eigen::VectorXd globalJointLowLimits) {

    int numSheaths = q.rows() / numConfigParams;
    // joints in this case means the configuration parameters of the system
    Eigen::VectorXd toReturn(numSheaths * numConfigParams); // holds secondary task that will be projected onto null space of Jacobian
    toReturn.setZero();

    Eigen::VectorXd weights(numSheaths * numConfigParams);
    weights << 1.0, 1.0, 10.0, 1.0, 1.0, 10.0, 1.0, 1.0, 10.0;

    // calculate the secondary task for each configuration param
    // the secondary task for this system is to keep the joint (config param) values
    // as far away from the joint limits as possible -> therefore, keep the joint values
    // close to the midpoint of the joint limits
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (((i + 1) % 3) == 2) { // calculation for phi; not currently controlled
            toReturn(i) = 0.0;
        }
        else { // calculation for theta and length
            toReturn(i) = weights(i) * (-1.0 / (1.0 * (numSheaths * numConfigParams))) *
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
/*Eigen::VectorXd calculateJointVelocities(Eigen::MatrixXd J,
    Eigen::VectorXd tipVelocityDesired, Eigen::VectorXd nullspaceObjectiveVector) {
    double toler = 0.002; // used in rank estimation for pinv

    int numSheaths = J.cols() / numConfigParams;
    // calculate pseudoinverse of Jacobian
    Eigen::MatrixXd pinvJ = pinv(J, toler);
    // calculate image space
    Eigen::VectorXd imageSpace = pinvJ * tipVelocityDesired;
    // calculate nullspace
    Eigen::VectorXd nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
        redundancyResolutionOn *
        nullspaceObjectiveVector;
    //Vec nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
    //    redundancyResolutionOn *
    //    nullspaceObjectiveVector;
    // calculate joint velocities
    Eigen::VectorXd qDot = imageSpace + nullSpace;

    bool skip = false;
    if ((nullSpace * nullSpace.transpose()).sum() == 0) {
        skip = true;
    }

    bool nullSpaceOverspeed = false;
    Eigen::VectorXd nullspaceScaling = Eigen::VectorXd::Zero((numSheaths * numConfigParams_)); // scaling factor for downscaling the nullspace

    if (!skip) {
        for (int i = 0; i < (numSheaths * numConfigParams_); ++i) { // check every joint to see if it exceeds the designated max speed
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) { // calculate nullspace scaling factor
                double nullSpaceElementIdealValue = sgn(qDot(i)) * maxJointSpeed(i) - imageSpace(i);
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
    Eigen::VectorXd qDotScaling = Eigen::VectorXd::Zero((numSheaths * numConfigParams));
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) {
            qDotOverspeed = true;
        }
    }

    if (qDotOverspeed) {
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) {
                qDotScaling(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        qDot = qDot / qDotScaling.maxCoeff(); // apply recalculated scaling factors
    }

    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000001) { // somehow still overspeeding
            std::cout << "overspeeding" << std::endl;
        }
    }

    return qDot;
}
*/
Eigen::VectorXd calculateJointVelocities(Eigen::MatrixXd J,
    Eigen::VectorXd tipVelocityDesired, Eigen::VectorXd nullspaceObjectiveVector, std::vector<double> lockedjoints, int exceededJointTracker) {
    double toler = 0.002; // used in rank estimation for pinv

    int numSheaths = J.cols() / numConfigParams;
    // calculate pseudoinverse of Jacobian
    Eigen::MatrixXd pinvJ = pinv(J, toler);
    // calculate image space
    Eigen::VectorXd imageSpace = pinvJ * tipVelocityDesired;
    // calculate nullspace
    Eigen::VectorXd nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
        redundancyResolutionOn *
        nullspaceObjectiveVector;
    // calculate joint velocities
    Eigen::VectorXd qDot = imageSpace + nullSpace;

    bool skip = false;
    //if ((nullSpace * nullSpace.transpose()).sum() == 0) {
    //    skip = true;
    //}

    bool nullSpaceOverspeed = false;
    Eigen::VectorXd nullspaceScaling = Eigen::VectorXd::Zero((numSheaths * numConfigParams_)); // scaling factor for downscaling the nullspace

    if (!skip) {
        for (int i = 0; i < (numSheaths * numConfigParams_); ++i) { // check every joint to see if it exceeds the designated max speed
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) { // calculate nullspace scaling factor
                double nullSpaceElementIdealValue = sgn(qDot(i)) * maxJointSpeed(i) - imageSpace(i);
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
    Eigen::VectorXd qDotScaling = Eigen::VectorXd::Zero((numSheaths * numConfigParams));
    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) {
            qDotOverspeed = true;
        }
    }

    if (qDotOverspeed) {
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001) {
                qDotScaling(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        qDot = qDot / qDotScaling.maxCoeff(); // apply recalculated scaling factors
    }

    for (int i = 0; i < (numSheaths * numConfigParams); ++i) {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000001) { // somehow still overspeeding
            std::cout << "overspeeding" << std::endl;
        }
    }

    return qDot;
}


Eigen::VectorXd inverseKinematics(Eigen::VectorXd& q, Eigen::VectorXd& qDot,
    Eigen::VectorXd &displacementCalc,
    Eigen::VectorXd tipVelocityDesired, int mode, int& outOfWorkspace) {
    // these must be hardcoded
    Eigen::Vector3d cushion; // deadband around joint limits
    cushion << 0.01, 0.0, 0.5; // order is theta, phi, length with units of rad, rad, mm

    int numSheaths = q.rows() / numConfigParams;
    // these must be hardcoded
    Eigen::VectorXd globalJointLowLimits(numConfigParams * numSheaths); // minimum allowable lengths/angles of the joints
    globalJointLowLimits << 0.0000000000000001, -1000.0 * M_PI, lengthLowLimit, -q(5) / minAllowedRadii, -1000.0 * M_PI, lengthLowLimit, -q(8) / minAllowedRadii, -1000.0 * M_PI, lengthLowLimit;

    // these must be hardcoded
    Eigen::VectorXd globalJointHighLimits(numConfigParams * numSheaths); // maximum allowable lengths/angles of the joints
    globalJointHighLimits << q(2) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit, q(5) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit, q(8) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit;
    
    // perform differential kinematics calculation to obtain Jacobian
    Eigen::MatrixXd J_normal = generateJacobian(q);

    Eigen::MatrixXd J_use; // Jacobian that will be used for robot control
    Eigen::VectorXd tipVelocityDesired_use; // desired input that will be used for robot control
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
            //locked_joints.push_back(7);
            J_use = J_normal.block(0, 0, 5, numConfigParams_ * numSheaths);
            tipVelocityDesired_use = tipVelocityDesired.head(5);
            //std::cout << "hello" << std::endl;
            //Eigen::MatrixXd temp1(4, numConfigParams_ * numSheaths);
            //temp1 << J_normal.row(0), J_normal.row(1), J_normal.row(2), J_normal.row(4);
            //J_use = temp1;

            //Eigen::VectorXd temp2(4);
            //temp2 << tipVelocityDesired(0), tipVelocityDesired(1), tipVelocityDesired(2), tipVelocityDesired(4);
            //tipVelocityDesired_use = temp2;
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
    Eigen::VectorXd nullspaceObjectiveVector = generateNullSpaceObjectiveVector(q,
        jointMidPoints,
        globalJointHighLimits,
        globalJointLowLimits);

    for (int i : locked_joints) {
        J_use.col(i) = Eigen::VectorXd::Zero(J_use.rows());
        nullspaceObjectiveVector(i) = 0.0;

        switch (i % numConfigParams) {
        case 1: // phi, units in rad
            if (i == 1) {
                globalJointLowLimits(i) = M_PI;
                globalJointHighLimits(i) = M_PI;
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

    Eigen::VectorXd qNew = q + qDot * deltaT; // recalculate config params

    int loop_counter = 0;
    Eigen::VectorXd joints_exceeded(numSheaths * numConfigParams);
    joints_exceeded.setZero();

    if (applyJointLimits) {
        bool finishedChecking = false;
        int numJoints = numSheaths * numConfigParams;
        Eigen::VectorXd limitHighJoint = globalJointHighLimits;
        Eigen::VectorXd limitLowJoint = globalJointLowLimits;

        while (!finishedChecking) {
            // calculate joint limits dynamically
            for (int i = 0; i < numSheaths; ++i) {
                if (i > 0) {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numConfigParams), -1.0 * (q(i * numConfigParams + 2) / minAllowedRadii)); // theta
                }
                else {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numConfigParams), 0.0); // theta
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
                J_use.col(exceededJoint) = Eigen::VectorXd::Zero(J_use.rows()); // wipe the column in the Jacobian of the current param
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

std::pair<Eigen::Vector3d, Eigen::Vector3d> calculate_circleTrajectory_position_velocity_f(Eigen::Vector3d start_point,
    double elapsed_time,
    double tip_speed,
    int out_of_workspace,
    Eigen::Vector3d pos_cur) {
    double radius = 20;
    Eigen::Vector3d euler_angles(0, M_PI / 4, 0);
    Eigen::Vector3d circle_center_point(start_point(0) - radius,
        start_point(1),
        start_point(2));

    //std::cout << "euler_angles: " << euler_angles << std::endl;
    //std::cout << "circle_center_point: " << circle_center_point << std::endl;

    // Calculate position 
    double omega = tip_speed / radius;
    double theta = omega * elapsed_time;

    //std::cout << "omega: " << omega << std::endl;
    //std::cout << "theta: " << theta << std::endl;

    double x = radius * std::cos(theta) + circle_center_point(0);
    double y = radius * std::sin(theta) + circle_center_point(1);
    double z = circle_center_point(2);

    //std::cout << "x: " << x << std::endl;
    //std::cout << "y: " << y << std::endl;
    //std::cout << "z: " << z << std::endl;

    // Calculate velocity
    double vx = -tip_speed * std::sin(theta);
    double vy = tip_speed * std::cos(theta);
    double vz = 0.0;

    // Rotation Matrices
    Eigen::Matrix3d Rz = getRotz_rad(euler_angles(2));

    Eigen::Matrix3d Ry = getRoty_rad(euler_angles(1));

    Eigen::Matrix3d Rx = getRotx_rad(euler_angles(0));

    //std::cout << "Rx: " << Rx << std::endl;
    //std::cout << "Ry: " << Ry << std::endl;
    //std::cout << "Rz: " << Rz << std::endl;

    // Final Rotation Matrix
    Eigen::Matrix3d rotation_matrix = Rz * Ry * Rx;

    // Rotate Velocity
    Eigen::Vector3d velocity = rotation_matrix * Eigen::Vector3d(vx, vy, vz);

    // Output position
    Eigen::Vector3d position = rotation_matrix * (Eigen::Vector3d(x, y, z) - start_point) + start_point;

    if (out_of_workspace == 1) {
        velocity(0) = (position(0) - pos_cur(0)) / deltaT;
        velocity(1) = (position(1) - pos_cur(1)) / deltaT;
        velocity(2) = (position(2) - pos_cur(2)) / deltaT;
    }

    // Return results using std::pair
    return std::make_pair(position, velocity);
}

double solve_equation(double f, double a, double b, double c) {
    // Initial guess (you may need to adjust this based on your specific problem)
    double th = -10.0;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double func = a + b * std::atan(th) + c * th - f;
        double derivative = b / (1 + th * th) + c;

        double delta = func / derivative;
        th -= delta;

        if (std::abs(delta) < EPSILON) {
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

Eigen::VectorXd createArray() {
    const int size = 500;
    Eigen::VectorXd arr(size);

    // Create a sequence from 1 to 500
    Eigen::VectorXd sequence = Eigen::VectorXd::LinSpaced(size, 1, size);

    // Calculate the sine part
    arr = 10 * (2 * M_PI * sequence / size).array().sin();

    // Generate random values and apply conditions
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (int i = 0; i < size; ++i) {
        double random_value1 = dis(gen);
        double random_value2 = dis(gen);

        if (random_value1 > 0.95) {
            arr(i) += 20;
        }
        if (random_value2 < 0.05) {
            arr(i) -= 20;
        }
    }

    return arr;
}

double sortTan(double prev, double input) {
    double a = abs(prev - input);
    double b = abs(prev - (input + M_PI));
    double c = abs(prev - (input - M_PI));
    double d = abs(prev - (input + 2 * M_PI));
    double e = abs(prev - (input - 2 * M_PI));

    double toReturn = a;
    double actualReturn = input;

    if (toReturn > b) {
        toReturn = b;
        actualReturn = input + M_PI;
    }

    if (toReturn > c) {
        toReturn = c;
        actualReturn = input - M_PI;
    }

    if (toReturn > d) {
        toReturn = d;
        actualReturn = input + 2 * M_PI;
    }

    if (toReturn > e) {
        toReturn = e;
        actualReturn = input - 2 * M_PI;
    }

    return actualReturn;
}

Eigen::Vector3d rotToEulXYZ(Eigen::Matrix3d rot, Eigen::Vector3d prev_ang) {
    Eigen::Vector3d toReturn = Eigen::Vector3d::Zero();

    double diff;

    double y1 = std::asin(rot(0, 2));
    double x1 = std::atan2(rot(1, 2) / (-std::cos(y1)), rot(2, 2) / std::cos(y1));
    double z1 = std::atan2(rot(0, 1) / (-std::cos(y1)), rot(0, 0) / std::cos(y1));
    //double diff1 = (std::pow((prev_ang(0) - x1), 2) + std::pow((prev_ang(1) - y1), 2) + std::pow((prev_ang(2) - z1), 2));
    diff = (std::pow(prev_ang(0) - sortTan(prev_ang(0), x1), 2) + std::pow(prev_ang(1)-y1, 2) + std::pow(prev_ang(2) - sortTan(prev_ang(2), z1), 2));
    toReturn << sortTan(prev_ang(0), x1), y1, sortTan(prev_ang(2), z1);

    double y2 = sgn(y1) * M_PI - y1;
    double x2 = std::atan2(rot(1, 2) / (-std::cos(y2)), rot(2, 2) / std::cos(y2));
    double z2 = std::atan2(rot(0, 1) / (-std::cos(y2)), rot(0, 0) / std::cos(y2));
    double diff2 = (std::pow(prev_ang(0) - sortTan(prev_ang(0), x2), 2) + std::pow(prev_ang(1) - y2, 2) + std::pow(prev_ang(2) - sortTan(prev_ang(2), z2), 2));
    if (diff2 < diff) {
        toReturn << sortTan(prev_ang(0), x2), y2, sortTan(prev_ang(2), z2);
        diff = diff2;
    }
   
    double y3 = sgn(y1) * -M_PI - y1;
    double x3 = std::atan2(rot(1, 2) / (-std::cos(y3)), rot(2, 2) / std::cos(y3));
    double z3 = std::atan2(rot(0, 1) / (-std::cos(y3)), rot(0, 0) / std::cos(y3));
    double diff3 = (std::pow(prev_ang(0) - sortTan(prev_ang(0), x3), 2) + std::pow(prev_ang(1) - y3, 2) + std::pow(prev_ang(2) - sortTan(prev_ang(2), z3), 2));
    if (diff3 < diff) {
        toReturn << sortTan(prev_ang(0), x3), y3, sortTan(prev_ang(2), z3);
        diff = diff3;
    }

    double y4 = sgn(y1) * -2* M_PI + y1;
    double x4 = std::atan2(rot(1, 2) / (-std::cos(y4)), rot(2, 2) / std::cos(y4));
    double z4 = std::atan2(rot(0, 1) / (-std::cos(y4)), rot(0, 0) / std::cos(y4));
    double diff4 = (std::pow(prev_ang(0) - sortTan(prev_ang(0), x4), 2) + std::pow(prev_ang(1) - y4, 2) + std::pow(prev_ang(2) - sortTan(prev_ang(2), z4), 2));
    if (diff4 < diff) {
        toReturn << sortTan(prev_ang(0), x4), y4, sortTan(prev_ang(2), z4);
    }

    return toReturn;
}

Eigen::Vector3d wrap2rel2PI(Eigen::Vector3d ang, Eigen::Vector3d init_ang) {
    Eigen::Vector3d toReturn; // Initialize the output vector
    for (int i = 0; i < 3; ++i) {
        toReturn(i) = (ang(i) > (init_ang(i) + M_PI)) ? (ang(i) - 2 * M_PI) : ((ang(i) < (init_ang(i) - M_PI)) ? (ang(i) + 2 * M_PI) : (ang(i)));
    }
    return toReturn;
}

Eigen::MatrixXd generateJacobian2(Eigen::VectorXd q) {

    int numSheaths = q.rows() / numConfigParams_;

    Eigen::VectorXd theta(numSheaths);
    Eigen::VectorXd phi(numSheaths);
    Eigen::VectorXd length(numSheaths);

    // separate and store config params
    for (int i = 0; i < numSheaths; ++i) {
        theta(i) = q(i * numConfigParams_);
        phi(i) = q(i * numConfigParams_ + 1);
        length(i) = q(i * numConfigParams_ + 2);
    }

    Eigen::MatrixXd J(6, numSheaths * numConfigParams_); // Jacobian matrix, 6x3n
    J.setZero();

    Eigen::Vector3d tip_pos; // position of tip of robot
    Eigen::Matrix3d temp_matrix = Eigen::Matrix3d::Identity(); // placeholder matrix
    forwardKinematics(q, tip_pos, temp_matrix); // get global robot tip position

    for (int i = 0; i < numSheaths; ++i) {

        Eigen::Vector3d prev_tip_pos = { 0.0, 0.0, 0.0 }; // tip position of previous segment chain
        Eigen::Matrix3d prev_rot = Eigen::Matrix3d::Identity(); // orientation matrix of tip of previous segment chain

        if (i > 0) {
            forwardKinematics(q.segment(0, i * numConfigParams_), prev_tip_pos, prev_rot); // saves tip position and orientation of previous segment
        }

        Eigen::Vector3d pos_ith_tip = { 0.0, 0.0, 0.0 }; // tip position of current segment
        Eigen::Matrix3d temp_matrix = Eigen::Matrix3d::Identity(); // Not used in this function call
        forwardKinematics(q.segment(0, (i + 1) * (numConfigParams_)), pos_ith_tip, temp_matrix); // saves tip position of current segment

        /** Jacobian formulation **/
        Eigen::MatrixXd M_i(6, 6);
        M_i << Eigen::Matrix3d::Identity(), -1 * skew(tip_pos - pos_ith_tip),
            Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

        Eigen::Vector3d J_v_i_correspondingToTheta = prev_rot * getRotz_rad(phi(i)) *
            Eigen::Vector3d(
                length(i) * (-1.0 * std::pow((theta(i) + (1e-6)), -2) + std::pow((theta(i) + (1e-6)), -2) * std::cos(theta(i)) + std::pow((theta(i) + (1e-6)), -1) * std::sin(theta(i))),
                0.0,
                length(i) * (-1.0 * std::pow((theta(i) + (1e-6)), -2) * std::sin(theta(i)) + std::pow((theta(i) + (1e-6)), -1) * std::cos(theta(i))));

        Eigen::Vector3d J_v_i_correspondingToPhi = prev_rot * (length(i) / (theta(i) + (1e-6))) * Eigen::Vector3d(-1.0 * std::sin(phi(i)) * (1.0 - std::cos(theta(i))), std::cos(phi(i)) * (1.0 - std::cos(theta(i))), 0.0);

        Eigen::Vector3d J_v_i_correspondingToLength = prev_rot * getRotz_rad(phi(i)) * Eigen::Vector3d(std::pow((theta(i) + (1e-6)), -1) * (1.0 - std::cos(theta(i))), 0.0, std::pow((theta(i) + (1e-6)), -1) * std::sin(theta(i)));

        Eigen::MatrixXd J_v_i(3, numConfigParams_); // holds all upper components of the Jacobian
        J_v_i << J_v_i_correspondingToTheta, J_v_i_correspondingToPhi, J_v_i_correspondingToLength;

        Eigen::Vector3d e2(0.0, 1.0, 0.0);
        Eigen::Vector3d e3(0.0, 0.0, 1.0);
        Eigen::MatrixXd J_w_i(3, numConfigParams_);
        J_w_i << prev_rot * getRotz_rad(phi(i)) * e2, prev_rot* e3, Eigen::Vector3d::Zero();

        Eigen::MatrixXd J_i(6, numConfigParams_);
        J_i << J_v_i, J_w_i;

        // performs adjoint transformation M_i (6x6) on
        // each sub-Jacobian J_i (6x3) and places the
        // resulting 6x3 matrix into the J matrix
        J.block(0, i * numConfigParams_, 6, numConfigParams_) = M_i * J_i;

    }

    return J;
}

void calculateJointVelocities2(Eigen::MatrixXd J,
    Eigen::VectorXd tipVelocityDesired, Eigen::VectorXd nullspaceObjectiveVector,
    Eigen::VectorXd &qDot, Eigen::VectorXd &imagespace, Eigen::VectorXd &nullspace) {
    double toler = 0.002; // used in rank estimation for pinv

    int numSheaths = J.cols() / numConfigParams_;
    // calculate pseudoinverse of Jacobian
    Eigen::MatrixXd pinvJ = J.transpose() *
        pinv(J * J.transpose() + (1e-6) * Eigen::MatrixXd::Identity(J.rows(), J.rows()), toler);
    // calculate image space
    Eigen::VectorXd imageSpace = pinvJ * tipVelocityDesired;
    // calculate nullspace
    Eigen::VectorXd nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams_), (numSheaths * numConfigParams_)) - pinvJ * J) *
        redundancyResolutionOn *
        nullspaceObjectiveVector;
    // calculate joint velocities
    qDot = imageSpace + nullSpace;
    imagespace = imageSpace;
    nullspace = nullSpace;
}

void readCSVMatrix(const std::string& filename, std::map<int, std::vector<double>>& toReturn) {
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
        std::vector<double> input;
        int count = 0;
        for (const auto& value : row) {
            input.push_back(value);
            count++;
        }
        std::cout << count << std::endl;
        toReturn[othercount] = input;
        othercount++;
    }
}

void readCSVMultiline(const std::string& filename, std::map<int, std::vector<double>>& toReturn) {
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
        std::vector<double> input;
        int count = 0;
        for (const auto& value : row) {
            input.push_back(value);
            count++;
        }
        toReturn[othercount] = input;
        othercount++;
    }
}

int main()
{
    //fesetround(FE_TONEAREST); // Round to nearest (even) like MATLAB

    /* theta_init_ << deg2rad(25), deg2rad(25), 0.00005; // initial bending angles of each sheath
    phi_init_ << 180.0 * M_PI / 180.0, -90.0 * M_PI / 180.0, 0.0 * M_PI / 180.0; // initial bending plane offsets of each sheath
    length_init_ << 10.0, 10.0, 10.0; // initial length of the steerable section of each sheath in mm
    delta_ << 2.3, 1.8, 1.0; // radius (in meters) of cross sectional circle centered at sheath centerline and passes through the tendon locations

    Eigen::VectorXd q(numConfigParams * numTubes);
    Eigen::VectorXd qDot(numConfigParams * numTubes);

    jointMidPoints << M_PI / 4.0, 0, 25.0, 0.0, 0.0, 25.0, 0.0, 0.0, 25.0;

    maxJointSpeed << M_PI / 6.0, 1000.0 * M_PI, 6.0,
        M_PI / 6.0, 1000.0 * M_PI, 6.0,
        M_PI / 6.0, 1000.0 * M_PI, 6.0;

    Eigen::MatrixXd tendonAngOffsets(numTotalTubes, numTendonsPerTube);
    tendonAngOffsets.row(0) << -1.0 * M_PI, -1.0 * 3.0 * M_PI / 2.0; // outer sheath
    tendonAngOffsets.row(1) << -1.0 * M_PI / 4.0, -1.0 * 3.0 * M_PI / 4.0; // intermediate sheath
    tendonAngOffsets.row(2) << -1.0 * 0.0, -1.0 * M_PI / 2.0; // TODO: inner sheath (currently unused)

    for (int i = 0; i < numTotalTubes; ++i) {
        tendonGeometry[i] << std::cos(tendonAngOffsets(i, 0)),
            std::cos(tendonAngOffsets(i, 1)),
            std::sin(tendonAngOffsets(i, 0)),
            std::sin(tendonAngOffsets(i, 1));

        invTendonGeometry[i] = tendonGeometry[i].inverse();
    }

    Eigen::VectorXd q2(9);
    Eigen::VectorXd xdot(3);
    Eigen::VectorXd ns(9);

    q2 << 0.206490867,	3.141592654,	23.872,	1.207449508, -1.615709922,	15.228,	5.00E-05,	0,	0.516;
    xdot << -0.795772196,	0.574894668,	0.364703211;
    ns << 0.001144692,	0, -0.003203685,	0, -0.003210294,	0.002324936,	0,	0,	0;

    Eigen::MatrixXd J_normal = generateJacobian2(q2);
    Eigen::MatrixXd J_use;
    J_use = J_normal.block(0, 0, 3, numConfigParams_ * 3);

    Eigen::VectorXd qDot2(9);
    Eigen::VectorXd imagespace(9);
    Eigen::VectorXd nullspace(9);
        
    calculateJointVelocities2(J_use, xdot, ns, qDot2, imagespace, nullspace);

    std::cout << "qdot: " << qDot2 << "\n" << std::endl;
    std::cout << "imagespace: " << imagespace << "\n" << std::endl;
    std::cout << "nullspace: " << nullspace << "\n" << std::endl;

    Eigen::MatrixXd A(2,3);

    A << 1e-6, 0, 0, 1e-12, 0, 0;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = svd.singularValues();

    if (!s.isZero(0)) {
        int r1 = (s.array() > 0.002).count(); // estimate effective rank
        r1 = std::max(r1, 1);
        std::cout << "r1: " << r1 << std::endl;

        Eigen::MatrixXd U = svd.matrixU();
        Eigen::MatrixXd V = svd.matrixV();

        std::cout << "U: " << U << "\n" << std::endl;
        std::cout << "V: " << V << "\n" << std::endl;
        std::cout << "s: " << s << "\n" << std::endl;

        // V.rightCols(V.cols() - r1).setZero();
        V.conservativeResize(V.rows(), r1);
        // U.rightCols(U.cols() - r1).setZero();
        U.conservativeResize(U.rows(), r1);
        s.conservativeResize(r1);

        s = s.array().inverse();
        std::cout << "s inverse: " << s << "\n" << std::endl;
        std::cout << "s diagonal: " << Eigen::MatrixXd(s.asDiagonal()) << "\n" << std::endl;

        //Eigen::VectorXd test(s.size());
        //for (int i = 0; i < test.size(); ++i) {
        //    if (s(i) > 0.8) {
        //        test(i) = 1 / s(i);
        //    }
        //    else {
        //        test(i) = 0;
        //    }
        //}

        std::cout << V * s.asDiagonal() * U.transpose() << std::endl;
        std::cout << A * (V * s.asDiagonal() * U.transpose()) << std::endl;
        //std::cout << V * test.asDiagonal() * U.transpose() << std::endl;
    }

    std::cout << q2.array().inverse() << std::endl; */

    //for (int i = 0; i < numTubes; ++i) {
    //    q_sim.segment((i * numConfigParams), numConfigParams) << theta_init_(i), phi_init_(i), length_init_(i);
    //}

    //std::vector<std::vector<double>> test;
    //readCSVDub("C:\\Users\\ch256744\\BCH Dropbox\\Phillip Tran\\ExtendedICRAPaper\\more experiments 1\\test2.csv", test);

    //std::vector<double> desired;
    //std::vector<double> actual;
    //std::vector<double> actual2;
    //for (int i = 0; i < test.size(); ++i) {
    //    for (int j = 0; j < test[i].size(); ++j) {
    //        //std::cout << test[i][j] << std::endl;
    //        if (j == 3) {
    //            desired.push_back(test[i][j]);
    //        }

    //        if (j == 9) {
    //            actual.push_back(test[i][j]);
    //        }

    //        if (j == 8) {
    //            actual2.push_back(test[i][j]);
    //        }
    //    }
    //}

    //MedianFilter mf = MedianFilter(11);
    //MedianFilter mf2 = MedianFilter(11);

    //std::vector<MedianFilter> medfilt;
    //for (int i = 0; i < 2; ++i) {
    //    medfilt.push_back(MedianFilter(11));
    //}

    //std::vector<SlidingModeFilter> slidingfilt;
    //for (int i = 0; i < 2; ++i) {
    //    slidingfilt.push_back(SlidingModeFilter(1, 1, 1, 0.01));
    //}

    //std::vector<double> temp;
    //std::vector<double> temp2;
    //std::vector<double> temp3;

    //for (int i = 0; i < desired.size(); ++i) {
    //   temp.push_back(mf.filter(desired[i]));
    //   //temp2.push_back(medfilt[0].filter(actual[i]));
    //   //temp3.push_back(medfilt[1].filter(actual2[i]));
    //   temp2.push_back(slidingfilt[0].filter(actual[i]));
    //   temp3.push_back(slidingfilt[1].filter(actual2[i]));
    //}

    ////writeCSV("C:/Users/ch256744/BCH Dropbox/Phillip Tran/ExtendedICRAPaper/AbdulsCode/MatlabCode/medfilttest63.csv", temp3);
    ////writeCSV("C:/Users/ch256744/BCH Dropbox/Phillip Tran/ExtendedICRAPaper/AbdulsCode/MatlabCode/medfilttest62.csv", temp2);

    //Eigen::Matrix3d dummyRotMM_i_wrtG1 = Eigen::Matrix3d::Identity();
    //forwardKinematics(q_sim, pTip_init_sim, dummyRotMM_i_wrtG1);

    ////for (int j = 0; j < 5000; ++j) {
    ////    auto time = j * deltaT;
    ////    Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity();
    ////    
    ////    forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);

    ////    std::pair<Eigen::Vector3d, Eigen::Vector3d> result =
    ////        calculate_circleTrajectory_position_velocity_f(pTip_init_sim, time, 0.5, outOfWorkspace, pTip_sim);
    ////    automated_velocity = result.second;

    ////    Eigen::VectorXd tipVelocityDesired(6); // command from manual control
    ////    tipVelocityDesired << automated_velocity(0), automated_velocity(1), automated_velocity(2), 0.0, 0.7 * ((2*M_PI) / 500.0) * std::cos(j * ((2 * M_PI) / 500.0)), 0.7 * ((2 * M_PI) / 500.0) * std::cos(j * ((2 * M_PI) / 500.0));

    ////    Eigen::VectorXd displacementCalc(nOfDrivenMotors);
    ////    displacementCalc = Eigen::VectorXd::Zero(nOfDrivenMotors); // motor displacement

    ////    //writeVectorToCSV("C:/Users/ch256744/Downloads/q_output_precise.csv", q_sim);
    ////    //writeVectorToCSV("C:/Users/ch256744/Downloads/qdot_output_precise.csv", qDot_sim);
    ////    //writeVectorToCSV("C:/Users/ch256744/Downloads/tipspeed_output_precise.csv", tipVelocityDesired);
    ////    
    ////    q_sim = inverseKinematics(q_sim, qDot_sim, displacement_sim, tipVelocityDesired, operatingMode, outOfWorkspace);
    ////    Eigen::MatrixXd jacobian = generateJacobian(q_sim);


    ////    forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);
    ////    writeVectorToCSV("C:/Users/ch256744/BCH Dropbox/Phillip Tran/ExtendedICRAPaper/AbdulsCode/MatlabCode/unlock3.csv", pTip_sim);
    ////    unemployed += outOfWorkspace;
    ////}

    //Eigen::MatrixXd t = Eigen::MatrixXd::Identity(4, 4);
    //std::cout << t << std::endl;

    //Eigen::VectorXd b(4);
    //b << 1, 2, 3, 4;
    //std::cout << b.head(3) << std::endl;
    //Eigen::MatrixXd n(4,4);
    //n << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    //std::cout << n.block(0,0,3,3).transpose() << std::endl;

    //auto p = getRotx(-177.34) * getRoty(-142.842) * getRotz(0.0024);
    //Eigen::Vector3d c;
    //c << deg2rad(2.64946), deg2rad(-37.1817), deg2rad(179.998);
    //Eigen::Vector3d u;
    //u << deg2rad(-2.023), deg2rad(-46.8271), deg2rad(177.134);
    ////auto g = wrap2rel2PI(rotToEulXYZ(p, c), u);
    //auto g = wrap2rel2PI(rotToEulXYZ(p, c), u);
    //std::cout << "Ans: " << rad2deg(g(0)) << " " << rad2deg(g(1)) << " " << rad2deg(g(2)) << std::endl;
    //std::cout << std::cos(179.9998) << std::endl;

    //std::map<int, std::vector<double>> test;
    //readCSVMatrix("C:\\Users\\ch256744\\BCH Dropbox\\Phillip Tran\\ExtendedICRAPaper\\ActualExperiments\\PathFollowingWithKinematicRedundancy\\M.csv", test);
    //auto a = test[2];
    //std::cout << a[23000] << std::endl;

    //std::map<int, std::vector<double>> calib_holder;
    //std::vector<double> l1;
    //std::vector<double> l2;
    //std::vector<double> tdOS;
    //std::vector<double> td45;
    //std::vector<double> td135;
    //std::vector<double> td225;
    //std::vector<double> td315;
    //readCSVMultiline("C:\\Users\\ch256744\\BCH Dropbox\\Phillip Tran\\ExtendedICRAPaper\\more experiments 1\\CalibrationInputOS.csv", calib_holder);
    //l2 = calib_holder[0];
    //td45 = calib_holder[1];
    //td135 = calib_holder[2];
    //td225 = calib_holder[3];
    //td315 = calib_holder[4];
    //std::cout << l2[9900] << std::endl;
    //for (int i = l2.size() - 101; i < l2.size(); ++i) {
    //    std::cout << "l2: " << l2[i] << std::endl;
    //}

    //Eigen::Vector3d length;
    //length << 1, 20, 0;
    //double phi_pos = 0;
    //Eigen::Vector3d phi;
    //phi << 0, -2, 0;
    //int i = 1;
    //Eigen::Vector2d f;
    //f << 1, 1;

    //double p = 0.865 + 0.0067 * length(1);
    //double q = 2.235 - 0.0067 * length(1);
    //double r = 7.465 + 0.1267 * length(1);

    ////double sol1 = (f(0) * ((0.003818 * length(i) + 0.02101))) / std::cos((phi_pos + phi(i - 1)) + (1e-8));
    ////double sol2 = (f(1) * ((0.003818 * length(i) + 0.02101))) / std::sin((phi_pos + phi(i - 1)) + (1e-8));

    //double sol1 = f(0) / (std::cos(phi_pos + phi(i - 1) + (1e-8)));
    //sol1 = 6;
    //sol1 = p * sol1 * (std::abs(sol1) < q) + ((r * sol1) - (sgn(sol1)) * (q * (r - p))) * (std::abs(sol1) >= q);
    ///*sol1 = deg2rad(sol1);*/

    //std::cout << "sol1: " << sol1 << std::endl;

    double a = 4.0;
    std::cout << 3 * (a > 3.0) << std::endl;
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
