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
    return Eigen::MatrixXd(toReturn.pseudoInverse());
}

double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

double rad2deg(double rad) {
    return rad * 180.0 / M_PI;
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

/* void directKinematicsV2_f(const Eigen::VectorXd& q,
    Eigen::Vector3d& g_end_tip_pos,
    Eigen::Matrix3d& rotMM_i_wrtG)
{
    int numberOfSegments = 1;

    Eigen::VectorXd theta(numberOfSegments);
    Eigen::VectorXd phi(numberOfSegments);
    Eigen::VectorXd lengthh(numberOfSegments);

    // store the config space params in separate vectors
    for (int i = 0; i < numberOfSegments; ++i)
    {
        theta(i) = q(i * 3);
        phi(i) = q(i * 3 + 1);
        lengthh(i) = q(i * 3 + 2);
    }

    // relative tip positions of every segment (sheath)
    Eigen::VectorXd p_i_rel_x(numberOfSegments);
    Eigen::VectorXd p_i_rel_y(numberOfSegments);
    Eigen::VectorXd p_i_rel_z(numberOfSegments);

    // tip positions of every segment in global ref frame
    Eigen::VectorXd g_tip_pos_x(numberOfSegments);
    Eigen::VectorXd g_tip_pos_y(numberOfSegments);
    Eigen::VectorXd g_tip_pos_z(numberOfSegments);

    // calculate tip positions of every segment in global ref frame
    for (int i = 0; i < numberOfSegments; ++i)
    {
        if (std::abs(theta(i)) < 0.0000000000000001)
        { // prevent divide by 0 error
            theta(i) = 0.0000000000000001;
        }

        Eigen::Matrix3d rotM_i_wrtG = Eigen::Matrix3d::Identity();

        for (int ii = 0; ii < i; ++ii)
        { // calculate orientation matrices of previous sheath(s)
            Eigen::Matrix3d zRotMatrixx;
            zRotMatrixx << std::cos(phi(ii)), -std::sin(phi(ii)), 0.0,
                std::sin(phi(ii)), std::cos(phi(ii)), 0.0,
                0.0, 0.0, 1.0;

            Eigen::Matrix3d yRotMatrix;
            yRotMatrix << std::cos(theta(ii)), 0.0, std::sin(theta(ii)),
                0.0, 1.0, 0.0,
                -std::sin(theta(ii)), 0.0, std::cos(theta(ii));

            rotM_i_wrtG = rotM_i_wrtG * zRotMatrixx * yRotMatrix;
        }

        Eigen::Matrix3d zRotMatrix;
        // calculate rotation about current z-axis by phi(i)
        zRotMatrix << std::cos(phi(i)), -std::sin(phi(i)), 0.0,
            std::sin(phi(i)), std::cos(phi(i)), 0.0,
            0.0, 0.0, 1.0;

        Eigen::Vector3d tip_pos_rel;
        tip_pos_rel << (lengthh(i) / theta(i)) * (1.0 - std::cos(theta(i))),
            0.0,
            (lengthh(i) / theta(i))* std::sin(theta(i));

        p_i_rel_x(i) = (zRotMatrix * tip_pos_rel).x();
        p_i_rel_y(i) = (zRotMatrix * tip_pos_rel).y();
        p_i_rel_z(i) = (zRotMatrix * tip_pos_rel).z();

        Eigen::Vector3d p_i_rell = rotM_i_wrtG * Eigen::Vector3d(p_i_rel_x(i),
            p_i_rel_y(i),
            p_i_rel_z(i));

        if (i == 0)
        { // if first segment
            g_tip_pos_x(i) = p_i_rell.x();
            g_tip_pos_y(i) = p_i_rell.y();
            g_tip_pos_z(i) = p_i_rell.z();
        }
        else
        { // account for previous segments
            g_tip_pos_x(i) = g_tip_pos_x(i - 1) + p_i_rell.x();
            g_tip_pos_y(i) = g_tip_pos_y(i - 1) + p_i_rell.y();
            g_tip_pos_z(i) = g_tip_pos_z(i - 1) + p_i_rell.z();
        }
    }

    // store global position of end tip
    g_end_tip_pos << g_tip_pos_x(numberOfSegments - 1),
        g_tip_pos_y(numberOfSegments - 1),
        g_tip_pos_z(numberOfSegments - 1);

    // rotation matrix of tip of robot wrt ref frame 0 (global)
    // rotMM_i_wrtG = Eigen::Matrix3d::Identity();

    for (int i = 0; i < numberOfSegments; ++i)
    { // calculate rotation transformation of robot tip wrt global frame
        Eigen::Matrix3d zRotMMatrixx;
        zRotMMatrixx << std::cos(phi(i)), -std::sin(phi(i)), 0.0,
            std::sin(phi(i)), std::cos(phi(i)), 0.0,
            0.0, 0.0, 1.0;

        Eigen::Matrix3d yRotMMatrix;
        yRotMMatrix << std::cos(theta(i)), 0.0, std::sin(theta(i)),
            0.0, 1.0, 0.0,
            -std::sin(theta(i)), 0.0, std::cos(theta(i));

        rotMM_i_wrtG = rotMM_i_wrtG * zRotMMatrixx * yRotMMatrix;
    }
} */

/* Eigen::MatrixXd differentialKinematicsV1_f(const Eigen::VectorXd& q)
{
    int numberOfSegments = 1;

    Eigen::VectorXd theta(numberOfSegments);
    Eigen::VectorXd phi(numberOfSegments);
    Eigen::VectorXd lengthh(numberOfSegments);

    // separate and store config params
    for (int i = 0; i < numberOfSegments; ++i)
    {
        theta(i) = q(i * 3);
        phi(i) = q(i * 3 + 1);
        lengthh(i) = q(i * 3 + 2);
    }

    Eigen::MatrixXd J(6, numberOfSegments * 3); // Jacobian matrix, 6x3n
    Eigen::Vector3d pTip;                       // position of tip of robot
    Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity(); // Not used in this function call
    directKinematicsV2_f(q, pTip, dummyRotMM_i_wrtG); // saves global robot tip position

    for (int i = 0; i < numberOfSegments; ++i)
    {
        Eigen::Vector3d posPre; // tip position of previous segment chain
        Eigen::Matrix3d rotPre = Eigen::Matrix3d::Identity(); // orientation matrix of tip of previous segment chain

        if (i == 0)
        { // no previous segment(s)
            posPre.setZero();
        }
        else
        {
            Eigen::VectorXd qq(3 * i);
            for (int ii = 0; ii < i; ++ii)
            { // stores config params of previous segment(s)
                qq.segment(ii * 3, 3) << theta(ii), phi(ii), lengthh(ii);
            }

            directKinematicsV2_f(qq, posPre, rotPre); // saves tip position and orientation of previous segment
        }

        Eigen::Matrix3d zRotMMatrix; // R_z(phi)
        zRotMMatrix << std::cos(phi(i)), -std::sin(phi(i)), 0.0,
            std::sin(phi(i)), std::cos(phi(i)), 0.0,
            0.0, 0.0, 1.0;

        Eigen::Vector3d e2(0.0, 1.0, 0.0);
        Eigen::Vector3d e3(0.0, 0.0, 1.0);

        Eigen::VectorXd qqq(3 * (i + 1));
        for (int iii = 0; iii <= i; ++iii)
        { // store config params of all segments up to and including current segment
            qqq.segment(iii * 3, 3) << theta(iii), phi(iii), lengthh(iii);
        }

        Eigen::Vector3d pos_ith_tip; // tip position of current segment
        Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity(); // Not used in this function call
        directKinematicsV2_f(qqq, pos_ith_tip, dummyRotMM_i_wrtG); // saves tip position of current segment

        Eigen::MatrixXd M_i(6, 6);
        M_i << Eigen::Matrix3d::Identity(), -1 * skew(pTip - pos_ith_tip),
            Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

        Eigen::Vector3d J_v_i_correspondingToTheta = rotPre * zRotMMatrix *
            Eigen::Vector3d(
                lengthh(i) * (-std::pow(theta(i), -2) + std::pow(theta(i), -2) * std::cos(theta(i)) + std::pow(theta(i), -1) * std::sin(theta(i))),
                0,
                lengthh(i) * (-std::pow(theta(i), -2) * std::sin(theta(i)) + std::pow(theta(i), -1) * std::cos(theta(i))));

        Eigen::Vector3d J_v_i_correspondingToPhi = rotPre * (lengthh(i) / theta(i)) * Eigen::Vector3d(
            -1.0 * std::sin(phi(i)) * (1.0 - std::cos(theta(i))),
            std::cos(phi(i)) * (1.0 - std::cos(theta(i))),
            0.0);

        Eigen::Vector3d J_v_i_correspondingToLengthh = rotPre * zRotMMatrix * std::pow(theta(i), -1) *
            Eigen::Vector3d(
                1.0 - std::cos(theta(i)),
                0.0,
                std::sin(theta(i)));

        Eigen::MatrixXd J_v_i(3, 3); // holds all upper components of the Jacobian
        J_v_i << J_v_i_correspondingToTheta, J_v_i_correspondingToPhi, J_v_i_correspondingToLengthh;

        Eigen::MatrixXd J_w_i(3, 3);
        J_w_i << rotPre * zRotMMatrix * e2, rotPre* e3, Eigen::Vector3d::Zero();

        Eigen::MatrixXd J_i(6, 3);
        J_i << J_v_i, J_w_i;

        J.block(0, i * 3, 6, 3) = M_i * J_i;
    }

    return J;
} */

/* Eigen::VectorXd qDotCalculate(const Eigen::MatrixXd& J,
    const Eigen::VectorXd& ksiDesired, int totalNumberOfJoints,
    double redunRes, const Eigen::VectorXd& omegaDifWRTq,
    const Eigen::VectorXd& maxJointSpeed)
{
    double toler = 0.002; // used in rank estimation for pinv

    // calculate pseudoinverse of Jacobian
    Eigen::MatrixXd pinvJ = pinv(J, toler);
    std::cout << "pinvj\n" << pinvJ << std::endl;
    // calculate image space
    Eigen::VectorXd imageSpace = pinvJ * ksiDesired;
    // calculate nullspace
    Eigen::VectorXd nullSpace = (Eigen::MatrixXd::Identity(totalNumberOfJoints, totalNumberOfJoints) - pinvJ * J) * redunRes * omegaDifWRTq;
    // calculate joint velocities
    Eigen::VectorXd qDot = imageSpace + nullSpace;

    

    bool overSpeed = false;
    Eigen::VectorXd overspeedScales = Eigen::VectorXd::Zero(totalNumberOfJoints); // scaling factor for downscaling the nullspace
    for (int i = 0; i < totalNumberOfJoints; ++i)
    { // check every joint to see if it exceeds the designated max speed
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
        { // calculate nullspace scaling factor
          // TODO: unclear why the ideal value for the secondary task is equal to the maximum joint speed - the image space
            double nullSpaceElementIdealValue = (qDot(i) / std::abs(qDot(i))) * maxJointSpeed(i) - imageSpace(i);
            overspeedScales(i) = nullSpace(i) / nullSpaceElementIdealValue;
            overSpeed = true;
        }
    }

    if (overSpeed)
    { // TODO: check if the logic makes sense for this scaling
      // TODO: what happens if overspeedScales only has negative values? that means that the below code doesn't run
        if (overspeedScales.maxCoeff() > std::numeric_limits<double>::epsilon())
        {                                                             // check that maximum scaling factor is greater than 0
            qDot = imageSpace + nullSpace / overspeedScales.maxCoeff(); // downscale nullspace
        }
    }

    // Check again if new qDot is beyond the speed limit
    bool overSpeedTwo = false;
    Eigen::VectorXd overspeedScalesTwo = Eigen::VectorXd::Zero(totalNumberOfJoints);
    for (int i = 0; i < totalNumberOfJoints; ++i)
    {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
        {
            overSpeedTwo = true;
        }
    }

    if (overSpeedTwo)
    { // TODO: check if the logic makes sense for this scaling
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < totalNumberOfJoints; ++i)
        {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
            {
                overspeedScalesTwo(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        // TODO: check if the logic makes sense for this scaling
        qDot = qDot / overspeedScalesTwo.maxCoeff(); // apply recalculated scaling factors
    }

    return qDot;
} */

/* Eigen::VectorXd calculateOmegaDifWRTq(int totalNumberOfJoints, const Eigen::VectorXd& q,
    const Eigen::VectorXd& jointMidPoints,
    const Eigen::VectorXd& jointHighLimits,
    const Eigen::VectorXd& jointLowLimits,
    double minAllowedRadii)
{
    Eigen::VectorXd omegaDifWRTq(totalNumberOfJoints);
    omegaDifWRTq.setZero();

    for (int ii = 0; ii < totalNumberOfJoints; ++ii)
    {
        if (((ii + 1) % 3) == 2)
        {
            omegaDifWRTq(ii) = 0.0;
        }
        else
        {
            omegaDifWRTq(ii) = (-1.0 / (2.0 * totalNumberOfJoints)) *
                (q(ii) - jointMidPoints(ii)) /
                ((jointHighLimits(ii) - jointLowLimits(ii)) *
                    (jointHighLimits(ii) - jointLowLimits(ii)));
        }
    }
    return omegaDifWRTq;
} */

/* std::vector<double> inverseK(double x_act, double y_act, double z_act, double OSL, double OSA) {
    Eigen::Vector4d pos_act(x_act, y_act, z_act, 1.0);
    double OSA_r = deg2rad(OSA);
    double x, y, z;
    if (std::abs(OSA) > 0) {
        Eigen::MatrixXd trans_OS(4,4);
        trans_OS << std::cos(OSA_r), 0.0, std::sin(OSA_r), (OSL / (OSA_r))* (1 - std::cos(OSA_r)), 0.0, 1.0, 0.0, 0.0, -1.0 * std::sin(OSA_r), 0.0, std::cos(OSA_r), (OSL / (OSA_r))* std::sin(OSA_r), 0.0, 0.0, 0.0, 1.0;
        Eigen::VectorXd pos = trans_OS.inverse() * pos_act;
        x = pos(0);
        y = pos(1);
        z = pos(2);
    }
    else {
        x = x_act;
        y = y_act;
        z = z_act;
    }

    double phi_rad = std::atan2(deg2rad(y), deg2rad(x));
    double th_rad = M_PI - 2 * std::atan(z / std::sqrt(x * x + y * y));
    double l = z * (th_rad / std::sin(th_rad));
    Eigen::Matrix2d beta;
    beta << std::cos(deg2rad(45)), std::cos(deg2rad(135)), std::sin(deg2rad(45)), std::sin(deg2rad(135));
    Eigen::Matrix2d beta_inv = beta.inverse();
    Eigen::Vector2d condenseVars((th_rad / 1.0) * std::cos(phi_rad) * 1.8 * 1.0, ((th_rad) / 1.0) * std::sin(phi_rad) * 1.8 * 1.0);
    Eigen::Vector2d ft = beta_inv * condenseVars;
    double td45 = ft(0);
    double td135 = ft(1);

    std::vector<double> toReturn;
    toReturn.push_back(td45);
    toReturn.push_back(td135);
    toReturn.push_back(rad2deg(phi_rad));
    toReturn.push_back(rad2deg(th_rad));
    toReturn.push_back(l);

    return toReturn;
} */

/* void positionToJointAngles(double px, double py, double pz, Eigen::VectorXd& toReturn) {

    auto p = std::sqrt(px * px + py * py + pz * pz);
    auto th2 = 2.0 * std::atan2(sqrt(py * py + px * px) / p, pz / p);
    auto phi2 = std::atan2(py / sqrt(py * py + px * px), px / sqrt(py * py + px * px));
    auto l = 0.0;

    if (phi2 <= 0) {
        phi2 = 2.0 * M_PI + phi2;
    }

    if (th2 == 0.0) {
        l = pz;
    }
    else {
        l = pz / (std::sin(th2) / th2);
    }

    auto th2_d = rad2deg(th2);
    auto phi2_d = rad2deg(phi2);

    int feasible = (l >= 0) & (l < 60) & (th2_d >= -120) & (th2_d <= 120);

    // if (feasible == 0) {
    // 	RCLCPP_INFO(this->get_logger(), "OUTER SHEATH REQUIRED --> th: %f and l: %f", th2_d, l);
    // 	RCLCPP_INFO(this->get_logger(), "\n");
    // }

    th2_d = std::min(120.0, std::max(-120.0, th2_d));
    l = std::min(60.0, std::max(0.0, l));

    toReturn << th2_d, phi2_d, l, 0.0, 0.0, feasible;
} */

/* Eigen::Matrix4d getTransformMat(Eigen::Matrix3d rotMat, double x, double y, double z) {
    Eigen::Matrix4d toReturn;
    toReturn << rotMat(0, 0), rotMat(0, 1), rotMat(0, 2), x,
        rotMat(1, 0), rotMat(1, 1), rotMat(1, 2), y,
        rotMat(2, 0), rotMat(2, 1), rotMat(2, 2), z,
        0.0, 0.0, 0.0, 1.0;

    return toReturn;
} */

/* std::vector<double> transformToJointVariables(Eigen::Matrix4d transform) {
    int feasible = 0;
    auto phi2 = 0.0;
    auto th2 = 0.0;
    auto th1 = 0.0;
    auto l2 = 0.0;
    auto l1 = 0.0;
    auto a = 0.0, b = 0.0, c = 0.0, d = 0.0;

    if (transform(1, 3) >= 0) {
        phi2 = std::acos(transform(1, 1));
    } else {
        phi2 = 2 * M_PI - std::acos(transform(1, 1));
    }

    std::cout << "phi2 " << phi2 << "\n" << std::endl;

    th2 = std::atan2(transform(1, 2) / std::sin(phi2), transform(1, 0) / std::sin(phi2));
    th1 = std::atan2(transform(2, 1) / std::sin(phi2), -1.0 * transform(0, 1) / std::sin(phi2));

    std::cout << "th2 " << th2 << "\n" << std::endl;
    std::cout << "th1 " << th1 << "\n" << std::endl;

    if (th2 != 0.0) {
        if (phi2 != 0.0) {
            l2 = (transform(1, 3) * (th2)) / ((1 - std::cos(th2)) * std::sin(phi2));

            if (th1 == 0) {
                l1 = (transform(2, 3) - 
                    (l2 / th2) * 
                    (std::sin(th2) * std::cos(th1) + sin(th1) * std::cos(phi2) * (std::cos(th2) - 1)));
            } else {
                l1 = (transform(2, 3) -
                    (l2 / th2) *
                    (std::sin(th2) * std::cos(th1) + sin(th1) * std::cos(phi2) * (std::cos(th2) - 1))) / 
                    (std::sin(th1) / th1);
            }
        } else {
            if (th1 == 0.0) {
                a = 0.0;
                b = (1 - std::cos(th2)) / th2;
                c = 1.0;
                d = std::sin(th2) / th2;
            } else {
                a = (1 - std::cos(th1)) / th1;
                b = (std::cos(th1) * (1 - std::cos(th2)) + std::sin(th1) * std::sin(th2)) / th2;
                c = std::sin(th1) / th1;
                d = (std::sin(th1) * (std::cos(th2) - 1) + std::cos(th1) * std::sin(th2)) / th2;
            }

            Eigen::Matrix2d toInv;
            toInv << a, b, c, d;
            Eigen::Vector2d toMult;
            toMult << transform(0, 3), transform(2, 3);
            auto result = toInv.inverse() * toMult;
            l1 = result(0);
            l2 = result(1);
        }
    } else {
        if (th1 == 0.0) {
            a = 0.0;
            c = 1.0;
        } else {
            a = (1 - std::cos(th1)) / th1;
            c = std::sin(th1) / th1;
        }
        b = std::sin(th1);
        d = std::cos(th1);
        Eigen::Matrix2d toInv;
        toInv << a, b, c, d;
        Eigen::Vector2d toMult;
        toMult << transform(0, 3), transform(2, 3);
        auto result = toInv.inverse() * toMult;
        l1 = result(0);
        l2 = result(1);
    }

    auto th1d = rad2deg(th1);
    auto th2d = rad2deg(th2);
    auto phi2d = rad2deg(phi2);

    std::vector<double> toReturn;
    toReturn.push_back(th1d);
    toReturn.push_back(l1);
    toReturn.push_back(th2d);
    toReturn.push_back(phi2d);
    toReturn.push_back(l2);
    toReturn.push_back(feasible);

    return toReturn;
} */

/* void decomposeTransform(Eigen::Matrix4d transform, Eigen::Matrix3d& rotation, Eigen::Vector3d& position) {

    rotation << transform(0, 0), transform(0, 1), transform(0, 2),
        transform(1, 0), transform(1, 1), transform(1, 2),
        transform(2, 0), transform(2, 1), transform(2, 2);

    position << transform(0, 3), transform(1, 3), transform(2, 3);
} */

/* void testOverwrite(Eigen::VectorXd& toReturn) {
    toReturn << 1, 2, 3, 4;
} */

/* double relativeAngleBetweenZAxes(Eigen::Vector3d angles, Eigen::Matrix3d rotCompare) {
    auto th_des = angles(0);
    auto phi_des = angles(1);

    auto rot_endtip = getRotz(phi_des) * getRoty(th_des);
    auto endtip_x = rot_endtip.col(0);
    auto endtip_z = rot_endtip.col(2);

    auto normal = endtip_x.cross(endtip_z);

    auto rotCompare_z = rotCompare.col(2);
    auto rotCompare_zproj = (rotCompare_z - normal * (rotCompare_z.dot(normal))).normalized();

    

    auto turn = rad2deg(std::acos(rotCompare_zproj.dot(endtip_z)));

    std::cout << turn << std::endl;
    //std::cout << rotCompare_zproj << std::endl;
    //std::cout << endtip_z << std::endl;

    auto cross_temp = endtip_z.cross(rotCompare_zproj);

    if (cross_temp.dot(normal) > 0.5) {
        turn *= -1;
    }

    return turn;
} */

/* void transformToJointVariables(Eigen::Matrix4d transform, Eigen::VectorXd& toReturn) {

    auto phi2 = 0.0;
    auto th2 = 0.0;
    auto th1 = 0.0;
    auto l2 = 0.0;
    auto l1 = 0.0;
    auto a = 0.0, b = 0.0, c = 0.0, d = 0.0;

    if (transform(1, 3) >= 0) {
        phi2 = std::acos(transform(1, 1));
    }
    else {
        phi2 = 2 * M_PI - std::acos(transform(1, 1));
    }

    th2 = std::atan2(transform(1, 2) / std::sin(phi2), transform(1, 0) / std::sin(phi2));
    th1 = std::atan2(transform(2, 1) / std::sin(phi2), -1.0 * transform(0, 1) / std::sin(phi2));

    if (th2 != 0.0) {
        if (phi2 != 0.0) {
            if (th1 == 0) {
                a = 0.0;
                b = (std::cos(phi2) * (1.0 - std::cos(th2))) / th2;
                c = 1.0;
                d = (std::sin(th2)) / th2;
            }
            else {
                a = (1.0 - std::cos(th1)) / th1;
                b = (std::cos(th1) * std::cos(phi2) * (1.0 - std::cos(th2)) + std::sin(th1) * std::sin(th2)) / th2;
                c = std::sin(th1) / th1;
                d = (std::sin(th1) * std::cos(phi2) * (std::cos(th2) - 1.0) + std::cos(th1) * std::sin(th2)) / th2;
            }
            Eigen::Matrix2d toInv;
            toInv << a, b, c, d;
            Eigen::Vector2d toMult;
            toMult << transform(0, 3), transform(2, 3);
            auto result = toInv.inverse() * toMult;
            l1 = result(0);
            l2 = result(1);
        }
        else {
            if (th1 == 0.0) {
                a = 0.0;
                b = (1.0 - std::cos(th2)) / th2;
                c = 1.0;
                d = std::sin(th2) / th2;
            }
            else {
                a = (1 - std::cos(th1)) / th1;
                b = (std::cos(th1) * (1.0 - std::cos(th2)) + std::sin(th1) * std::sin(th2)) / th2;
                c = std::sin(th1) / th1;
                d = (std::sin(th1) * (std::cos(th2) - 1.0) + std::cos(th1) * std::sin(th2)) / th2;
            }
            Eigen::Matrix2d toInv;
            toInv << a, b, c, d;
            Eigen::Vector2d toMult;
            toMult << transform(0, 3), transform(2, 3);
            auto result = toInv.inverse() * toMult;
            l1 = result(0);
            l2 = result(1);
        }
    }
    else {
        if (th1 == 0.0) {
            a = 0.0;
            c = 1.0;
        }
        else {
            a = (1.0 - std::cos(th1)) / th1;
            c = std::sin(th1) / th1;
        }
        b = std::sin(th1);
        d = std::cos(th1);
        Eigen::Matrix2d toInv;
        toInv << a, b, c, d;
        Eigen::Vector2d toMult;
        toMult << transform(0, 3), transform(2, 3);
        auto result = toInv.inverse() * toMult;
        l1 = result(0);
        l2 = result(1);
    }

    auto th1d = rad2deg(th1);
    auto th2d = rad2deg(th2);
    auto phi2d = rad2deg(phi2);

    int feasible = (l1 >= 0) & (th1 > -90) & (th1 <= 0) &
        (l2 >= 0) & (l2 < 50) & (th2 >= -90) & (th2 <= 90);

    toReturn << th2d, phi2d, l2, th1d, l1, feasible;
} */

/* Eigen::Matrix4d buildTransform(Eigen::Matrix3d rotation, Eigen::VectorXd position) {

    Eigen::Matrix4d toReturn;
    toReturn << rotation(0, 0), rotation(0, 1), rotation(0, 2), position(0),
        rotation(1, 0), rotation(1, 1), rotation(1, 2), position(1),
        rotation(2, 0), rotation(2, 1), rotation(2, 2), position(2),
        0.0, 0.0, 0.0, 1.0;

    return toReturn;
} */

/* void forwardKinematics(const Eigen::VectorXd q,
    Eigen::Vector3d& end_tip_pos,
    Eigen::Matrix3d& rotationMatrix)
{
    int numSheaths = q.rows() / numConfigParams;

    Eigen::VectorXd theta(numSheaths);
    Eigen::VectorXd phi(numSheaths);
    Eigen::VectorXd length(numSheaths);

    // store the config space params in separate vectors
    for (int i = 0; i < numSheaths; ++i)
    {
        theta(i) = (q(i * numConfigParams) == 0.0) ? std::numeric_limits<double>::epsilon() : q(i * numConfigParams);
        phi(i) = q(i * numConfigParams + 1);
        length(i) = q(i * numConfigParams + 2);
    }

    rotationMatrix = Eigen::Matrix3d::Identity();
    end_tip_pos.setZero();

    for (int k = 0; k < numSheaths; ++k) {
        rotationMatrix = Eigen::Matrix3d::Identity();

        for (int j = 0; j < (k + 1); ++j) {

            if (j < k) {
                rotationMatrix = rotationMatrix * (getRotz_rad(phi(j))) * (getRoty_rad(theta(j)));
            }
            else {
                Eigen::Vector3d tip_pos_rel;
                tip_pos_rel << (length(j) / theta(j)) * (1.0 - std::cos(theta(j))),
                    0.0,
                    (length(j) / theta(j))* std::sin(theta(j));
                end_tip_pos = end_tip_pos + rotationMatrix * (getRotz_rad(phi(j))) * tip_pos_rel;
            }

            if (j == (numSheaths - 1)) {
                rotationMatrix = rotationMatrix * (getRotz_rad(phi(j))) * (getRoty_rad(theta(j)));
            }
        }
    }
} */

/* Eigen::MatrixXd generateJacobian(Eigen::VectorXd q) {
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
    Eigen::Vector3d pTip;  // position of tip of robot
    // this is required b/c in some uses of directKinematicsV2_f
    // the orientation matrix is saved but not in this instance
    Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity(); // Not used in this function call
    forwardKinematics(q, pTip, dummyRotMM_i_wrtG); // saves global robot tip position

    for (int i = 0; i < numSheaths; ++i) {

        Eigen::Vector3d posPre = { 0.0, 0.0, 0.0 }; // tip position of previous segment chain
        Eigen::Matrix3d rotPre = Eigen::Matrix3d::Identity(); // orientation matrix of tip of previous segment chain

        if (i > 0) {
            forwardKinematics(q.segment(0, i * numSheaths), posPre, rotPre); // saves tip position and orientation of previous segment
        }
        //std::cout << "posepre: " << posPre << std::endl;
        Eigen::Vector3d pos_ith_tip; // tip position of current segment
        // this is required b/c in some uses of directKinematicsV2_f
        // the orientation matrix is saved but not in this instance
        Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity(); // Not used in this function call
        forwardKinematics(q.segment(0, (i + 1) * (numSheaths)), pos_ith_tip, dummyRotMM_i_wrtG); // saves tip position of current segment
        //std::cout << "posithtip: " << pos_ith_tip << std::endl;

        Eigen::MatrixXd M_i(6, 6);
        M_i << Eigen::Matrix3d::Identity(), -1 * skew(pTip - pos_ith_tip),
            Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

        Eigen::Vector3d J_v_i_correspondingToTheta = rotPre * getRotz_rad(phi(i)) *
            Eigen::Vector3d(
                length(i) * (-1.0 * std::pow(theta(i), -2) + std::pow(theta(i), -2) * std::cos(theta(i)) + std::pow(theta(i), -1) * std::sin(theta(i))),
                0.0,
                length(i) * (-1.0 * std::pow(theta(i), -2) * std::sin(theta(i)) + std::pow(theta(i), -1) * std::cos(theta(i))));

        Eigen::Vector3d J_v_i_correspondingToPhi = rotPre * (length(i) / theta(i)) * Eigen::Vector3d(
            -1.0 * std::sin(phi(i)) * (1.0 - std::cos(theta(i))),
            std::cos(phi(i)) * (1.0 - std::cos(theta(i))),
            0.0);

        Eigen::Vector3d J_v_i_correspondingToLength = rotPre * getRotz_rad(phi(i)) * Eigen::Vector3d(
            std::pow(theta(i), -1) * (1.0 - std::cos(theta(i))),
            0.0,
            std::pow(theta(i), -1) * std::sin(theta(i)));

        Eigen::MatrixXd J_v_i(3, numConfigParams); // holds all upper components of the Jacobian
        J_v_i << J_v_i_correspondingToTheta, J_v_i_correspondingToPhi, J_v_i_correspondingToLength;

        Eigen::Vector3d e2(0.0, 1.0, 0.0);
        Eigen::Vector3d e3(0.0, 0.0, 1.0);
        Eigen::MatrixXd J_w_i(3, numConfigParams);
        J_w_i << rotPre * getRotz_rad(phi(i)) * e2, rotPre* e3, Eigen::Vector3d::Zero();

        Eigen::MatrixXd J_i(6, numConfigParams);
        J_i << J_v_i, J_w_i;

        // performs adjoint transformation M_i (6x6) on
        // each sub-Jacobian J_i (6x3) and places the
        // resulting 6x3 matrix into the J matrix
        J.block(0, i * numConfigParams, 6, numConfigParams) = M_i * J_i;
    }

    return J;
} */

/* Eigen::VectorXd generateNullSpaceObjectiveVector(Eigen::VectorXd q,
    Eigen::VectorXd jointMidPoints,
    Eigen::VectorXd jointHighLimits,
    Eigen::VectorXd jointLowLimits)
{
    // joints in this case means the configuration parameters of the system
    Eigen::VectorXd toReturn(numTubes * numConfigParams); // holds secondary task that will be projected onto null space of Jacobian
    toReturn.setZero();

    Eigen::VectorXd weights(numTubes * numConfigParams);
    weights << 1.0, 1.0, 10.0, 1.0, 1.0, 10.0, 1.0, 1.0, 10.0;

    // calculate the secondary task for each configuration param
    // the secondary task for this system is to keep the joint (config param) values
    // as far away from the joint limits as possible -> therefore, keep the joint values
    // close to the midpoint of the joint limits
    for (int i = 0; i < (numTubes * numConfigParams); ++i) {
        if (((i + 1) % 3) == 2) { // calculation for phi; not currently controlled
            toReturn(i) = 0.0;
        }
        else { // calculation for theta and length
            toReturn(i) = weights(i) * (-1.0 / (1.0 * (numTubes * numConfigParams))) *
                (q(i) - jointMidPoints(i)) /
                ((jointHighLimits(i) - jointLowLimits(i)) *
                    (jointHighLimits(i) - jointLowLimits(i)));
        }
    }

    return toReturn;
} */

/* Eigen::VectorXd calculateJointVelocities(Eigen::MatrixXd J,
    Eigen::VectorXd tipVelocityDesired, Eigen::VectorXd nullspaceObjectiveVector) {
    double toler = 0.002; // used in rank estimation for pinv

    // calculate pseudoinverse of Jacobian
    Eigen::MatrixXd pinvJ = pinv(J, toler);
    
    // calculate image space
    Eigen::VectorXd imageSpace = pinvJ * tipVelocityDesired;
    std::cout << "inverse matrix: " << pinvJ << std::endl;
    // calculate nullspace
    Eigen::VectorXd nullSpace = (Eigen::MatrixXd::Identity((numTubes * numConfigParams), (numTubes * numConfigParams)) - pinvJ * J) * 1.0 * nullspaceObjectiveVector;
    // calculate joint velocities
    Eigen::VectorXd qDot = imageSpace + nullSpace;

    bool nullSpaceOverspeed = false;
    Eigen::VectorXd nullspaceScaling = Eigen::VectorXd::Zero((numTubes * numConfigParams)); // scaling factor for downscaling the nullspace


    Eigen::Vector3d maxJointSpeedSub(M_PI / 6.0, 1000.0 * M_PI, 0.003);
    Eigen::MatrixXd maxJointSpeed(numTubes * numConfigParams, 1);
    maxJointSpeed << maxJointSpeedSub, maxJointSpeedSub, maxJointSpeedSub;

    for (int i = 0; i < (numTubes * numConfigParams); ++i)
    { // check every joint to see if it exceeds the designated max speed
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
        { // calculate nullspace scaling factor
            double nullSpaceElementIdealValue = (qDot(i) / std::abs(qDot(i))) * maxJointSpeed(i) - imageSpace(i);
            nullspaceScaling(i) = nullSpace(i) / nullSpaceElementIdealValue;
            nullSpaceOverspeed = true;
        }
    }

    if (nullSpaceOverspeed)
    {
        if (nullspaceScaling.maxCoeff() > std::numeric_limits<double>::epsilon())
        {                                                                // check that maximum scaling factor is greater than 0
            qDot = imageSpace + nullSpace / nullspaceScaling.maxCoeff(); // downscale nullspace
        }
    }

    // Check again if new qDot is beyond the speed limit
    bool qDotOverspeed = false;
    Eigen::VectorXd qDotScaling = Eigen::VectorXd::Zero((numTubes * numConfigParams));
    for (int i = 0; i < (numTubes * numConfigParams); ++i)
    {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
        {
            qDotOverspeed = true;
        }
    }

    if (qDotOverspeed)
    {
        qDot = imageSpace + nullSpace; // calculate joint velocities using original image space and nullspace
        for (int i = 0; i < (numTubes * numConfigParams); ++i)
        {
            if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000000000001)
            {
                qDotScaling(i) = std::abs(qDot(i)) / maxJointSpeed(i);
            }
        }
        qDot = qDot / qDotScaling.maxCoeff(); // apply recalculated scaling factors
    }

    for (int i = 0; i < (numTubes * numConfigParams); ++i)
    {
        if (std::abs(qDot(i)) > maxJointSpeed(i) + 0.0000000001)
        { // somehow still overspeeding
            std::cout << "overspeeding" << std::endl;
        }
    }

    return qDot;
} */

/* void jointSpaceToConfigurationSpace(int numberOfSegments,
    const Eigen::VectorXd tendonTotalDisplacement,
    const Eigen::VectorXd lengthhBackbone,
    Eigen::VectorXd phi,
    Eigen::VectorXd theta) {
    // initialize displacement due to other segments
    Eigen::MatrixXd accumulationsForSections(numTendonsPerTube, numberOfSegments);
    accumulationsForSections.setZero();

    Eigen::MatrixXd curvaturePair(2, numberOfSegments); // curvature components K1 and K2
    Eigen::VectorXd curvature(numberOfSegments);        // curvature K where K = sqrt(K1^2 + K2^2)
    Eigen::VectorXd phiPre(numberOfSegments);           // previous bending plane angle

    curvaturePair.setZero();
    curvature.setZero();
    phiPre.setZero();

    Eigen::VectorXd delta_(numTubes);
    delta_ << 0.0023, 0.0018, 0.0010;

    // initialize angular positions of tendons wrt axial cross sectional x axis
    Eigen::MatrixXd beta(numTotalTubes, numTendonsPerTube);
    // this must be hardcoded
    beta.row(0) << M_PI, 3.0 * M_PI / 2.0;                // outer sheath
    beta.row(1) << M_PI / 4.0, 3.0 * M_PI / 4.0;          // intermediate sheath
    beta.row(2) << 0.0, M_PI / 2.0;                       // TODO: inner sheath (currently unused)

    // this must be hardcoded
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> capitalPiMatrixi_[3]; 
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> invCapitalPiMatrixi_[3];
    for (int i = 0; i < numTotalTubes; ++i)
    {
        capitalPiMatrixi_[i] << std::cos(beta(i, 0)),
            std::cos(beta(i, 1)),
            std::sin(beta(i, 0)),
            std::sin(beta(i, 1));

        invCapitalPiMatrixi_[i] = capitalPiMatrixi_[i].inverse();
    }

    Eigen::VectorXd theta_init_(numTubes);
    theta_init_ << 0.0000000000000001, 0.0000000000000001, 0.0000000000000001; // initial bending angles of each section, constant
    Eigen::VectorXd phi_init_(numTubes);
    phi_init_ << 180.0 * M_PI / 180.0, -90.0 * M_PI / 180.0, 0.0 * M_PI / 180.0;           // initial bending plane offsets of each section, constant

    for (int i = 0; i < numberOfSegments; ++i)
    {
        // v.segment<2>(i) grabs a subsection of vector v that contains 2 indices
        // inclusive of index i ie v.segment<2>(0) returns [v(0), v(1)]
        Eigen::VectorXd tendonIsolatedForThisSection =
            (tendonTotalDisplacement.segment(i * 4, numTendonsPerTube) - tendonTotalDisplacement.segment(i * 4 + 2, numTendonsPerTube) / 2 - accumulationsForSections.col(i));

        curvaturePair.col(i) = (1 / (delta_(i) * lengthhBackbone(i))) * capitalPiMatrixi_[i] * tendonIsolatedForThisSection;

        // propagate tendon displacement forward to obtain
        // accumulated effects for the next segment
        for (int ii = i + 1; ii < numberOfSegments; ++ii)
        {
            // invCapitalPiMatrixi_ is a global variable that is initialized in the constructor
            accumulationsForSections.col(ii) += invCapitalPiMatrixi_[ii] * delta_(ii) * lengthhBackbone(i) * curvaturePair.col(i);
        }

        // calculate curvature from curvature components
        curvature(i) = curvaturePair.col(i).norm();

        // adjust curvature values so that they will never be negative zero because that messes up atan2
        auto k0 = (curvaturePair(0, i) == 0.0) ? abs(curvaturePair(0, i)) : curvaturePair(0, i);
        auto k1 = (curvaturePair(1, i) == 0.0) ? abs(curvaturePair(1, i)) : curvaturePair(1, i);

        // enforce deadband of -0.1 to 0.1 for curvature components
        // TODO: this makes the signed zero checking in the above two lines completely useless but let's not worry about that for now
        // k0 = (std::abs(k0) < deadbandThreshold) ? 0.0 : k0;
        // k1 = (std::abs(k1) < deadbandThreshold) ? 0.0 : k1;

        // calculate bending plane angle from curvature components
        phiPre(i) = std::atan2(k1, k0);

        // calculate bending angle
        // theta(i) = theta_init_(i) + curvature(i) * lengthhBackbone(i);
        theta(i) = curvature(i) * lengthhBackbone(i);

        if (i == 0)
        {
            // TODO: This is not a mistake for the current implementation of the robot.
            // TODO: The robot cannot control phi(0) because it does not have a rotation DOF
            // TODO: for the outer tube. Therefore, phi(0) will always be equal to the
            // TODO: hardcoded initial value for phi(0). In a different configuration of
            // TODO: of the robot that can control phi(0), the following code should be used:

            phi(0) = phi_init_(0);
            std::cout << "phi" << i << ": " << rad2deg(phiPre(i)) << std::endl;

            // TODO: commenting this out to see if it changes anything
            // TODO: the logic behind this as it was explained does not make sense
            if (std::abs(phiPre(0) - phi_init_(0)) > 0.0001)
            {
                theta(0) = 0.0000000000000001;
            }
        }
        else
        {
            // calculate relative bending plane angle
            //          phi(i) = phi_init_(i) + 0*(phiPre(i) - phiPre(i - 1));
            phi(i) = phiPre(i) - phiPre(i - 1);

            std::cout << "phi" << i << ": " << rad2deg(phi(i)) << std::endl;
        }

        if (std::abs(theta(i)) < 0.0000000000000001)
        { // hardcode theta to a very small value to resolve jittering of theta values
            theta(i) = 0.0000000000000001;
        }
    }

    // TODO: these are hardcoded because currently
    // TODO: the inner sheath does not have tendons and therefore
    // TODO: is not actuated so the values are substituted in
    phi(2) = phi_init_(2);
    theta(2) = theta_init_(2);
} */

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
Eigen::VectorXd calculateJointVelocities(Eigen::MatrixXd J,
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
    //Vec nullSpace = (Eigen::MatrixXd::Identity((numSheaths * numConfigParams), (numSheaths * numConfigParams)) - pinvJ * J) *
    //    redundancyResolutionOn *
    //    nullspaceObjectiveVector;
    // calculate joint velocities
    Eigen::VectorXd qDot = imageSpace + nullSpace;

    bool skip = false;
    if ((nullSpace * nullSpace.transpose()).sum() == 0) {
        //for (int j = 0; j < nullSpace.size(); ++j) {
        //    nullSpace(j) = EPSILON;
        //}

        //for (int k : lockedjoints) {
        //    nullSpace(k) = 0.0;
        //}

        //if (exceededJointTracker != -1) {
        //    nullSpace(exceededJointTracker) = 0.0;
        //}
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


Eigen::VectorXd inverseKinematics(Eigen::VectorXd &q, Eigen::VectorXd &qDot,
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
            locked_joints.push_back(7);
            // J_use = J_normal.block(0, 0, 5, numConfigParams_ * numSheaths);
            // tipVelocityDesired_use = tipVelocityDesired.head(5);
            Eigen::MatrixXd temp1(4, numConfigParams_ * numSheaths);
            temp1 << J_normal.row(0), J_normal.row(1), J_normal.row(2), J_normal.row(4);
            J_use = temp1;

            Eigen::VectorXd temp2(4);
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

int main()
{
    /* std::cout << "Hello World!\n";

    Eigen::VectorXd q(9);
    q.segment(0, 3) << 4, 5, 6;
    std::cout << "This is the vector:\n" << q << std::endl;

    Eigen::Matrix<double, 3, 2> m;
    std::cout << "This is the size of the matrix:\n" << m.size() << std::endl;

    Eigen::MatrixXd acc(2, 3);
    acc <<  1, 1,
            2, 2,
            3, 3;
    std::cout << "This is acc:\n" << acc << std::endl;
    std::cout << "This is the column:\n" << acc.col(0) << std::endl;

    Eigen::VectorXd b(12);
    b << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    std::cout << "This is b:\n" << b << std::endl;
    std::cout << "This is size of b: " << b.size() << std::endl;

    std::cout << "This is segment 1:\n" << b.segment<2>(0) << std::endl;
    std::cout << "This is segment 2:\n" << b.segment<2>(2) << std::endl;

    std::cout << "This is head of lettuce:\n" << b.head(3) << std::endl;

    Eigen::VectorXd c = (b.segment<2>(0) - b.segment<2>(2)) / 2;
    std::cout << "This is c:\n" << c << std::endl;

    std::cout << "This is a number: " << q.rows() << std::endl;

    Eigen::Matrix3d d;
    std::cout << "This is size of matrix: " << d.size() << std::endl;

    Eigen::Matrix3d zRotMatrix;
    zRotMatrix << 1,    1,     0,
        1,              1,                   0,
        0,                  0,                  1;
    std::cout << "This is z: " << zRotMatrix<< std::endl;

    Eigen::Vector3d u;
    Eigen::Vector3d y;

    y << 1, 2, 3;
    u = y.head(3);

    std::cout << "This is u:\n" << u << std::endl;

    std::cout << "This is math: " << 2 % 3 << std::endl;

    Eigen::VectorXd p(3);
    p(0) = 4;
    p(1) = 8;
    std::cout << "This is another head of lettuce:\n" << p.head(3) << std::endl;

    Eigen::VectorXd l(12);
    Eigen::Vector2d n(1, 2);

    for (int i = 0; i < 3; ++i) {
        l.segment(4 * i, 2) = n;
    }

    std::cout << "This is a vector:\n" << l << std::endl;

    Eigen::Matrix<double, 14, 19> huge;

    std::cout << "This is huge:\n" << huge << std::endl;

    huge << Eigen::MatrixXd::Identity(14, 14),
        Eigen::MatrixXd::Zero(14, 5);

    std::cout << "This is huge2:\n" << huge << std::endl;

    std::vector<double> big = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0,
            0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    
    int count = 0;
    for (int i = 0; i < 14; i++) {
        for (int j = 0; j < 19; j++) {
            huge(i, j) = big.at(count);
            count++;
        }
    }

    std::cout << "This is huge3:\n" << huge << std::endl;


    std::bitset<16> bitset1{ 1 };
    std::cout << bitset1 << std::endl;

    std::istringstream converter("221C");
    uint16_t val;
    converter >> std::hex >> val;
    std::cout << "Hex to int is: " << val << std::endl;

    int testarray1[7] = { 20 };
    std::fill(testarray1, testarray1 + 7, 20);
    for (int j = 0; j < 7; j++) {
        std::cout << testarray1[j] << std::endl;
    }

    double vel = 5.168413;
    double dir = -1.0;
    int con = 10;
    int32_t theint = vel * dir * con;
    std::cout << "This is a number: " << theint << std::endl;

    int num1 = 9;
    std::cout << "Rounded number: " << std::ceil(num1 / 2.0) << std::endl; */

    /* Eigen::MatrixXd beta(3, 2);
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> capitalPiMatrixi_[3];
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> invCapitalPiMatrixi_[3];
    beta.row(0) << M_PI, 3 * M_PI / 2;          
    beta.row(1) << M_PI / 4, 3 * M_PI / 4;      
    beta.row(2) << 0, M_PI / 2;  
    for (int i = 0; i < 3; ++i)
    {
        capitalPiMatrixi_[i] << std::cos(beta(i, 0)),
            std::cos(beta(i, 1)),
            std::sin(beta(i, 0)),
            std::sin(beta(i, 1));

        invCapitalPiMatrixi_[i] = capitalPiMatrixi_[i].inverse();
    }
    std::cout << "This is inv cap pi:\n" << invCapitalPiMatrixi_[1] << std::endl; */

    /* Eigen::Matrix3d zRotMatrixx;
    zRotMatrixx << std::cos(M_PI / 4), -std::sin(M_PI / 4), 0,
        std::sin(M_PI / 4), std::cos(M_PI / 4), 0,
        0, 0, 1;

    std::cout << "This is the matrix:\n" << zRotMatrixx << std::endl; */

    /*Eigen::MatrixXd test(2, 3);
    test.col(0) << 3, 4;
    std::cout << "This is the matrix:\n" << -test << std::endl;

    Eigen::VectorXd q(9);
    q.setZero();
    std::cout << "Number of rows:\n" << q << std::endl;
    q.segment(3, 3) << 3, 4, 5;
    std::cout << "Number of rows:\n" << q << std::endl; */

    /* Eigen::Vector3d t;
    t << -4* (-std::pow(3, -2) + std::pow(3, -2) * std::cos(3) + std::pow(3, -1) * std::sin(3)), 
        0, 
        4* (-std::pow(3, -2) * std::sin(3) + std::pow(3, -1) * std::cos(3));
    std::cout << "Number of rows:\n" << t << std::endl;
    std::cout << "max:\n" << t.head(2) << std::endl;

    Eigen::Vector3d t1;
    double sca = 4 * std::pow(3, -2);
    t1 << sca*(std::cos(3) - 1 + 3 * std::sin(3)),
        0,
        sca * (-1*std::sin(3) + 3 * std::cos(3));
    std::cout << "Number of rows:\n" << t1 << std::endl;

    Eigen::MatrixXd h(6, 3);
    h << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity();
    std::cout << "Number:\n" << h << std::endl;

    Eigen::MatrixXd z(6,6);
    Eigen::Matrix3d y;
    y << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    z << Eigen::Matrix3d::Identity(), -1 * y,
        Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();;

    Eigen::MatrixXd J(6, 3* 3);
    J.block(0, 6, 6, 3) = z * h;
    std::cout << "Number:\n" << J << std::endl;
        
    Eigen::MatrixXd j(6,9);
    for (int i = 0; i < 6; ++i)
    {
        for (int k = 0; k < 9; ++k)
        {
            j(i, k) = i + k;
        }
    }
    std::cout << "Number:\n" << j << std::endl;
    std::cout << "inv matrix:\n" << pinv(j, 0.002) << std::endl;
    // std::cout << "other method:\n" << j.completeOrthogonalDecomposition().pseudoInverse() << std::endl;
    auto b = j.completeOrthogonalDecomposition();
    b.setThreshold(0.002);
    auto c = b.pseudoInverse();
    // std::cout << "another one:\n" << c << std::endl;

    Eigen::MatrixXd l(b.pseudoInverse());
    std::cout << "third one:\n" << l << std::endl;

    Eigen::MatrixXd g = Eigen::MatrixXd::Identity(9, 9) - Eigen::MatrixXd(b.pseudoInverse()) * j;
    Eigen::MatrixXd v = Eigen::MatrixXd::Identity(9, 9) - pinv(j, 0.002) * j;
    std::cout << "aha:\n" << g << std::endl;
    std::cout << "aha2:\n" << v << std::endl;

    double r = -3 / 10.0;
    std::cout << "hello: " << signbit(r) << std::endl; */

    /* std::cout << "hi" << std::endl;

    Eigen::VectorXd q2(9);
    q2 << 1.00000000000000e-16, 3.14159265358979, 0.0300000000000000, 1.00000000000000e-16, 0.0, 0.0300000000000000, 1.00000000000000e-16, 0.0, 0.0300000000000000;
    std::cout << q2 << std::endl;

    Eigen::MatrixXd m1 = differentialKinematicsV1_f(q2);

    std::cout << m1 << std::endl;

    std::cout << pinv(m1,0.002) << std::endl;

    Eigen::Vector3d g = { 180.0 * M_PI / 180.0, -90.0 * M_PI / 180.0, 0.0 * M_PI / 180.0 };
    std::cout << g << std::endl; */
    
    /* double testarray[3] = {1, 2 ,3};
    Eigen::Vector3d testvector = { 4, 5, 6 };
    Eigen::VectorXd testvector2(3);

    testvector2(0) = testarray[0] + testvector(0);
    // std::cout << testvector2(0) << std::endl;

    Eigen::Matrix<double, 2, 3> h;
    h << 1, 2,
        3, 4,
        5, 6;

    // std::cout << h << std::endl;
    // std::cout << h.col(1).y() << std::endl; */

    /* Eigen::MatrixXd beta(3, 2);
    beta.row(0) << M_PI, 3.0 * M_PI / 2.0;       
    beta.row(1) << M_PI / 4.0, 3.0 * M_PI / 4.0;  
    beta.row(2) << 0.0, M_PI / 2.0;                  

    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> capitalPiMatrixi_[3];               
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> invCapitalPiMatrixi_[3];

    for (int i = 0; i < 3; ++i)
    {
        capitalPiMatrixi_[i] << std::cos(beta(i, 0)),
            std::cos(beta(i, 1)),
            std::sin(beta(i, 0)),
            std::sin(beta(i, 1));

        invCapitalPiMatrixi_[i] = capitalPiMatrixi_[i].inverse();
    }

    // std::cout << capitalPiMatrixi_[1] << std::endl;
    // std::cout << invCapitalPiMatrixi_[1] << std::endl;

    Eigen::VectorXd tendonTotalDisplacement(12);
    Eigen::MatrixXd curvaturePair(2, 3); // curvature components K1 and K2
    Eigen::VectorXd curvature(3);           // curvature K where K = sqrt(K1^2 + K2^2)
    Eigen::VectorXd phiPre(3);              // previous bending plane angle

    tendonTotalDisplacement << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ;
    Eigen::VectorXd tendonIsolatedForThisSection =
        (tendonTotalDisplacement.segment<2>(0 * 4) - tendonTotalDisplacement.segment<2>(0 * 4 + 2)) / 2;
    curvaturePair.col(0) = (1 / (0.0023 * 0.03)) * capitalPiMatrixi_[0] * tendonIsolatedForThisSection;
    curvature(0) = curvaturePair.col(0).norm();
    phiPre(0) = atan2(curvaturePair(1, 0), curvaturePair(0, 0));

    // std::cout << tendonIsolatedForThisSection << std::endl;
    std::cout << curvaturePair.col(0) << std::endl;
    std::cout << "curve pair 0: " << curvaturePair(0, 0) << std::endl;
    std::cout << "curve pair 1: " << curvaturePair(1, 0) << std::endl;
    std::cout << "atan2 " << atan2(0, -0) << std::endl;
    std::cout << "wtf " << atan2(curvaturePair(1, 0), curvaturePair(0, 0)) << std::endl;
    auto a = curvaturePair(1, 0);
    auto b = curvaturePair(0, 0);
    std::cout << "this is a: " << a << std::endl;
    std::cout << "this is b: " << b << std::endl;
    std::cout << "this is atan2 ab: " << atan2(abs(a), abs(b)) << std::endl;
    std::cout << phiPre(0) << std::endl;

    auto c = (a == 0.0) ? abs(b) : b;
    std::cout << "this is c: " << c << std::endl;
    std::cout << "this is b: " << b << std::endl;
    std::cout << "this is atan2 ab: " << atan2(a, c) << std::endl;

    std::cout << "this is cosine: " << cos(b) << std::endl;
    std::cout << "this is sine: " << sin(b) << std::endl; */

    /* double minAllowedRadii = 0.032;
    double lengthLowLimit = 0.001;
    double lengthHighLimit = 0.050;
    Eigen::VectorXd jointMidPoints(9);
    jointMidPoints << M_PI / 4.0, 0, 0.025, 0.0, 0.0, 0.025, 0.0, 0.0, 0.025;

    Eigen::VectorXd q(9);
    q << 1.00000000000000e-16,	3.14159265358979, - 0.000840000000000001,	1.40852527676149, - 1.57079632679490,	0.0450600000000000,	1.00000000000000e-16,	0,	0.0524760000000000;

    Eigen::VectorXd jointLowLimits(9);
    jointLowLimits << 0.0000000000000001, -1000.0 * M_PI, lengthLowLimit, -q(5) / minAllowedRadii, -1000.0 * M_PI, lengthLowLimit, -q(8) / minAllowedRadii, -1000.0 * M_PI, lengthLowLimit;

    Eigen::VectorXd jointHighLimits(9);
    jointHighLimits << q(2) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit, q(5) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit, q(8) / minAllowedRadii, 1000.0 * M_PI, lengthHighLimit;

    Eigen::MatrixXd jacobianmatrix(3, 9);
    jacobianmatrix.row(0) << -0.0393287095693388, 0, 0, 8.19324999572998e-18, 0, -3.37234888855183e-17, 0, 0, 4.40575259891950e-17;
    jacobianmatrix.row(1) << 4.81637782886867e-18,	0,	0,	0.0211267923438870,	0,	0.594711649641283,	0,	0,	0.986542083973620;
    jacobianmatrix.row(2) << -8.82710983233251e-18,	0,	1, - 0.0689960706472038,	0,	0.701390747929602,	0,	0,	0.163507542789277;

    Eigen::MatrixXd totestmatrix = jacobianmatrix;
    

    Eigen::VectorXd omegaDifWRTq = calculateOmegaDifWRTq(9,
        q,
        jointMidPoints,
        jointHighLimits,
        jointLowLimits,
        minAllowedRadii);

    for (int ii : {1, 4, 6, 7})
    {
        jacobianmatrix.col(ii) = Eigen::VectorXd::Zero(jacobianmatrix.rows());
        omegaDifWRTq(ii) = 0.0;
    }

    Eigen::VectorXd ksiDesired(3);
    ksiDesired << 0.0, 0.0005, 0.0;

    Eigen::Vector3d maxJointSpeedSub(M_PI / 6.0, 1000.0 * M_PI, 0.003);
    Eigen::MatrixXd maxJointSpeed(9, 1);
    maxJointSpeed << maxJointSpeedSub, maxJointSpeedSub, maxJointSpeedSub;

    Eigen::VectorXd qDot = qDotCalculate(
        jacobianmatrix,
        ksiDesired,
        9,
        0.0,
        omegaDifWRTq,
        maxJointSpeed);

    // std::cout << qDot << std::endl;

    std::cout << "before\n" << jacobianmatrix << std::endl;

    jacobianmatrix.col(3) = Eigen::VectorXd::Zero(jacobianmatrix.rows());
    jacobianmatrix.col(0) = Eigen::VectorXd::Zero(jacobianmatrix.rows());
    jacobianmatrix.col(2) = Eigen::VectorXd::Zero(jacobianmatrix.rows());
    jacobianmatrix.col(5) = Eigen::VectorXd::Zero(jacobianmatrix.rows());
    jacobianmatrix.col(8) = Eigen::VectorXd::Zero(jacobianmatrix.rows());

    std::cout << "after\n" << jacobianmatrix << std::endl;

    qDot = qDotCalculate(
        jacobianmatrix,
        ksiDesired,
        9,
        0.0,
        omegaDifWRTq,
        maxJointSpeed);

    std::cout << qDot << std::endl;

    totestmatrix = totestmatrix * 0;
    std::cout << pinv(totestmatrix, 0.002) << std::endl;

    std::cout << "\n" << std::ceil(0.01) << std::endl;

    double b = 4.6;
    double c = 5.0;

    int d = b + c;

    std::cout << d << std::endl; */

    /* Eigen::VectorXd b(9);
    b << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::VectorXd g = b.segment(1, 8);
    auto n = g.sum();
    std::cout << g << std::endl;
    std::cout << n << std::endl;
    Eigen::Vector3d posPre = { 1.0, 0.0, 0.0 };
    std::cout << posPre << std::endl;
    posPre.setZero();
    std::cout << posPre << std::endl;

    double a = 2.89308776595231886830;
    printf("ans = %.50f\n", std::cos(a));

    std::cout << -std::pow(0.25, -2) << std::endl;

    Eigen::Matrix3d zRotMMatrix; // R_z(phi)
    zRotMMatrix << 1, 0, 0.0,
        0, 1, 0.0,
        0.0, 0.0, 1.0;

    Eigen::Vector3d test = zRotMMatrix * Eigen::Vector3d(0, 2, 0);
    std::cout << test << std::endl;

    Eigen::MatrixXd beta(3, 2);
    beta.row(0) << M_PI, 3.0 * M_PI / 2.0;                // outer sheath
    beta.row(1) << M_PI / 4.0, 3.0 * M_PI / 4.0;          // intermediate sheath
    beta.row(2) << 0.0, M_PI / 2.0;

    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> capitalPiMatrixi_[3];
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> mat[3];
    for (int i = 0; i < 3; ++i)
    {
        capitalPiMatrixi_[i] << std::cos(beta(i, 0)), 
            std::cos(beta(i, 1)), 
            std::sin(beta(i, 0)), 
            std::sin(beta(i, 1));

        mat[i] = capitalPiMatrixi_[i].inverse();
    }
    std::cout << mat[0] << std::endl;
    std::cout << mat[1] << std::endl;
    std::cout << mat[2] << std::endl;

    Eigen::VectorXd lowLimitHolder(9);
    lowLimitHolder << 1, 2, 3, 4,5,6,7,8,9;
    Eigen::VectorXd highLimitHolder(9);
    highLimitHolder << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Eigen::VectorXd jointLowLimits(9); // minimum allowable lengths/angles of the joints
    Eigen::VectorXd jointHighLimits(9); // maximum allowable lengths/angles of the joints

    for (int i = 0; i < 9; ++i)
    {
        jointLowLimits[i] = lowLimitHolder[i];
    }

    std::cout << jointLowLimits << std::endl;

    Eigen::VectorXd j(3) = { 1,2,3 }; */

    /* Eigen::VectorXd q(numTotalTubes* numConfigParams);
    q << 0.5640, -0.0000, 0.2997, 0.3089, 0.0619, 0.0234, 0.0001, -0.0590, 0.0002;
    Eigen::VectorXd ksiDesired(6);
    ksiDesired << 0, 0, 0.0005, 0, 0, 0;
    Eigen::MatrixXd J_normal = differentialKinematicsV1_f(q);
    std::cout << "j normal " << J_normal << std::endl;
    Eigen::MatrixXd J_use;          // Jacobian that will be used for robot control
    Eigen::VectorXd ksiDesired_use; // desired input that will be used for robot control
    J_use = J_normal.block(0, 0, 3, numConfigParams * numTubes);
    ksiDesired_use = ksiDesired.head(3);
    std::cout << "j use " << J_use << std::endl;
    std::cout << "ksi " << ksiDesired_use << std::endl;

    Eigen::VectorXd jointLowLimits(numConfigParams*numTubes); // minimum allowable lengths/angles of the joints
    jointLowLimits << 0.0000000000000001, -1000.0 * M_PI, 0.01;
    Eigen::VectorXd jointHighLimits(numConfigParams* numTubes); // maximum allowable lengths/angles of the joints
    jointHighLimits << q(2) / 0.032, 1000.0 * M_PI, 0.05;
    Eigen::VectorXd jointMidPoints(numTubes*numConfigParams);
    jointMidPoints << M_PI / 4.0, 0, 0.025;
    Eigen::VectorXd omegaDifWRTq = calculateOmegaDifWRTq(numTubes*numConfigParams,
        q,
        jointMidPoints,
        jointHighLimits,
        jointLowLimits,
        0.032);

    std::cout << "omegadif " << omegaDifWRTq << std::endl;

    Eigen::Vector3d maxJointSpeedSub(M_PI / 6.0, 1000.0 * M_PI, 0.003);
    Eigen::MatrixXd maxJointSpeed(numConfigParams*numTubes, 1);
    maxJointSpeed << maxJointSpeedSub;

    Eigen::VectorXd qDot = qDotCalculate(J_use, ksiDesired_use, numConfigParams*numTubes, 0.0, omegaDifWRTq, maxJointSpeed);

    std::cout << "qdot " << qDot << std::endl;

    Eigen::VectorXd t(3);
    t(0) = 3;
    t(1) = 4;
    t(2) = 9;
    std::cout << "t " << t << std::endl; */

    /*Eigen::VectorXd q(numTotalTubes* numConfigParams);
    q << 0.5640, -0.0000, 0.2997, 0.3089, 0.0619, 0.0234, 0.0001, -0.0590, 0.0002;

    std::cout << q.segment(0, numTotalTubes * numConfigParams) << std::endl;

    Eigen::Vector3d maxJointSpeedSub(M_PI / 6.0, 1000.0 * M_PI, 0.003);
    Eigen::VectorXd maxJointSpeed(numTubes* numConfigParams);
    maxJointSpeed = Eigen::VectorXd::Zero(numTubes * numConfigParams);

    for (int i = 0; i < numTubes; ++i)
    {
        maxJointSpeed << maxJointSpeedSub;
    }

    std::cout << maxJointSpeed << std::endl;

    Eigen::MatrixXd maxJointSpeed_sim(numTubes * numConfigParams, 1);
    maxJointSpeed_sim << maxJointSpeedSub, maxJointSpeedSub, maxJointSpeedSub;

    std::cout << maxJointSpeed_sim << std::endl; */

    /* Eigen::VectorXd tendonTotalDisplacement_(numTendonsPerTube* nOfTendonMotors);
    tendonTotalDisplacement_ << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    Eigen::MatrixXd accumulationsForSections(numTendonsPerTube, numTubes);
    accumulationsForSections.setZero();
    std::cout << accumulationsForSections << std::endl;

    Eigen::MatrixXd tendonIsolatedForThisSection(numTendonsPerTube, numTubes);
    for (int i = 0; i < numTubes; ++i)
    {
        tendonIsolatedForThisSection.col(i) =
            (tendonTotalDisplacement_.segment(i * 4,numTendonsPerTube) - tendonTotalDisplacement_.segment(i * 4 + 2, numTendonsPerTube)) / 2;

        std::cout << tendonTotalDisplacement_.segment(i * 4, numTendonsPerTube) << std::endl;
        std::cout << tendonTotalDisplacement_.segment(i * 4 + 2, numTendonsPerTube) << std::endl;
    }
    std::cout << tendonIsolatedForThisSection << std::endl;


    Eigen::MatrixXd test(numTendonsPerTube, numTubes);
    test.col(0) << 4, 5;
    std::cout << test << std::endl;

    std::cout << tendonTotalDisplacement_.size() << std::endl; */

    /* std::vector<double> test;
    test = inverseK(5.0, 2.0, 10.0, 2.0, -10.0);
    for (int i = 0; i < test.size(); ++i) {
        std::cout << test[i] << std::endl;
    }
    Eigen::VectorXd inverse_kinematics_results(5);
    inverse_kinematics_results << 1, 2, 3, 4, 5;
    // std::cout << inverse_kinematics_results << std::endl;
    inverse_kinematics_results << 5, 4, 3, 2, 1;
    std::cout << inverse_kinematics_results << std::endl;

    std::cout << 37.7e-3 << std::endl;

    std::cout << std::max(10.7793,0.0001) << std::endl; */

    // Eigen::Matrix3d test;
    // test << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    // std::cout << test << std::endl;

    /* std::vector<double> configs = positionToJointAngles(-3.0, -1.0, 3.0);
    //for (int i = 0; i < configs.size(); ++i) {
    //    std::cout << configs[i] << std::endl;
    //}
    //std::cout << "\n" << std::endl;
    auto baoa = 45.0;
    auto eaoa = 1.0;
    auto roty1 = getRoty(baoa*eaoa);
    auto rotz1 = getRotz(configs[1]);
    auto roty2 = getRoty(configs[0]);
    auto roty3 = getRoty(eaoa);
    auto X = -3.0;
    auto Y = -1.0;
    auto Z = 3.0;
    Eigen::Matrix4d transform;
    Eigen::Matrix3d finalRot = roty1 * rotz1 * roty2 * roty3;
    //std::cout << finalRot << std::endl;
    //std::cout << "\n" << std::endl;
    transform << finalRot(0,0), finalRot(0,1), finalRot(0,2), X,
        finalRot(1,0), finalRot(1,1), finalRot(1,2), Y,
        finalRot(2,0), finalRot(2,1), finalRot(2,2), Z,
        0.0, 0.0, 0.0, 1.0;
    //std::cout << transform << std::endl;
    std::vector<double> results1 = transformToJointVariables(transform);
    for (int i = 0; i < results1.size(); ++i) {
        std::cout << results1[i] << std::endl;
    }

    std::cout << "hi1" << std::endl;
    Eigen::VectorXd j(4);
    j << 8, 9, 10, 11;
    std::cout << j << std::endl;
    testOverwrite(j);
    std::cout << "hi3" << std::endl;
    std::cout << j << std::endl;

    Eigen::Matrix3d m = Eigen::Matrix3d::Identity();
    Eigen::Vector3d p;
    p << 1, 2, 3;
    decomposeTransform(transform, m, p);
    std::cout << m << std::endl;
    std::cout << p << std::endl;

    for (int i = 0; i < 9; ++i) {
        std::cout << m(i) << std::endl;
    }

    std::vector<double> g;
    for (int i = 0; i < 16; ++i) {
        g.push_back(transform(i));
    }
    Eigen::Matrix4d transformation;
    transformation << g[0], g[4], g[8], g[12],
        g[1], g[5], g[9], g[13],
        g[2], g[6], g[10], g[14],
        g[3], g[7], g[11], g[15];
    std::cout << transform << std::endl;
    std::cout << transformation << std::endl;
    std::cout << g.size() << std::endl;

    Eigen::Matrix3d o;
    o << 1, 0, 0, 0, -1, 0, 0, 0, -1;
    std::cout << o << std::endl;
    //for (int i = 0; i < 16; ++i) {
    //    std::cout << g[i] << std::endl;
    //} */

    /* Eigen::Vector3d angles = positionToJointAngles(1, 1, 40);
    
    Eigen::Matrix3d rotation;
    rotation << -0.9613, - 0.2049, 0.1839, -0.1265,  0.9220, 0.3660, -0.2446, 0.3286, -0.9123;

    std::cout << rad2deg(relativeAngleBetweenZAxes(angles, rotation)) << std::endl; */

    /*
    //Eigen::VectorXd g = Eigen::VectorXd::Zero(2);

    //Eigen::VectorXd b = Eigen::VectorXd::Zero(2);
    //b << 2, 1;
    //Eigen::VectorXd c = Eigen::VectorXd::Zero(2);
    //c << 9, 4;
    //g = b + c;
    ////std::cout << g << std::endl;

    //int feasible = 0;
    //auto l = -5;
    //auto th = 20;

    //feasible = (l >= 0) & (l < 50) & (th>=-90) & (th<=90);

    //std::cout << feasible << std::endl;


    std::map<std::pair<int, int>, Eigen::Matrix3d> map;

    //map[std::pair<int, int>(1, 0)] = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d s;
    s << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    //std::cout << s << std::endl;
    //map[std::pair<int, int>(2, 1)] = s;

    //std::cout << map[std::pair<int, int>(2, 1)] << std::endl;

    std::map<int, Eigen::Matrix3d> map_use;
    std::map<int, Eigen::Matrix3d> map_yaw;

    readCSVDub("C:\\Users\\ch256744\\BCH Dropbox\\Phillip Tran\\ICRA 2025\\test_csv2.csv", map_use);
    readCSVDub("C:\\Users\\ch256744\\BCH Dropbox\\Phillip Tran\\ICRA 2025\\test_csv3.csv", map_yaw);

    //int othercount = 0;
    // Access data like a 2D vector
    //for (const auto& row : data) {
    //    Eigen::Matrix3d input = Eigen::Matrix3d::Identity();
    //    int count = 0;
    //    for (const auto& value : row) {
    //        input(count) = value;
    //        count++;
    //        //std::cout << value << " ";
    //    }
    //    map[std::pair<int, int>(othercount, 0)] = input.transpose();
    //    //std::cout << std::endl;
    //    othercount++;
    //}

    //std::cout << map[std::pair<int, int>(2, 0)] << std::endl;

    auto px = -8.0;
    auto py = 18.0;
    auto pz = 36.0;

    Eigen::Matrix3d useme;
    Eigen::Matrix3d useyaw;
    Eigen::Vector4d pos_useme;
    pos_useme << px, py, pz, 1.0;
    useme = map_use[18];
    useyaw = map_yaw[5];
    std::cout << useme << std::endl;
    std::cout << useyaw << std::endl;

    auto transform = buildTransform(useme, pos_useme);
    Eigen::VectorXd dist_sheath_angles_end;
    dist_sheath_angles_end = Eigen::VectorXd::Zero(6);
    transformToJointVariables(transform, dist_sheath_angles_end);
    std::cout << "Sheath info -> TH2_d: " << dist_sheath_angles_end(0) <<
        "PHI2_d: " << dist_sheath_angles_end(1) << 
        "L2: " << dist_sheath_angles_end(2) << 
        "OSA: " << dist_sheath_angles_end(3) << 
        "OSL: " << dist_sheath_angles_end(4) << std::endl;
    */

    /* // use gimbal angles in calculations if corresponding button is held on stylus
    //dual_sheath_mode = (msg.effort[1]) ? true : false;
    //RCLCPP_INFO(this->get_logger(), "GIMBAL MODE STATUS: %d",
    //    int(dual_sheath_mode));
    //RCLCPP_INFO(this->get_logger(), "\n");

    //distal_pullback = (msg.effort[3]) ? true : false;

    // MIGHT NOT NEED THIS
    // if ((!msg.effort[2] && prox_advanced = false) || (!msg.effort[2] && prox_advanced = true)) {
    // 	run_continuous = true;
    // } else {
    // 	run_continuous = false;
    // }

    // UNCOMMENT FOR BACKUP PLAN
    //run_continuous = (!msg.effort[2]) ? true : false;

    // get position vector
    Eigen::VectorXd position(4);
    position << 1.0, 2.0, 3.0, 1.0;

    // get gimbal angles
    Eigen::Vector3d gimbals;
    gimbals << 0.707, 1.0, 0.0;

    auto px = position[0];
    auto py = position[1];
    auto pz = position[2];

    px = std::min(10.0, std::max(-40.0, px));
    py = std::min(25.0, std::max(-10.0, py));
    pz = std::min(40.0, std::max(1.0, pz));

    // convert gimbals to degrees
    auto rx = rad2deg(gimbals[0]);
    auto ry = rad2deg(gimbals[1]);
    auto rz = rad2deg(gimbals[2]);

    //RCLCPP_INFO(this->get_logger(), "RECEIVED POSITION----->  X:%f Y %f Z: %f", px, py, pz);
    //RCLCPP_INFO(this->get_logger(), "\n");
    //RCLCPP_INFO(this->get_logger(), "RECEIVED GIMBALS----->  X:%f Y %f Z: %f", rx, ry, rz);
    //RCLCPP_INFO(this->get_logger(), "\n");

    // checkLimits(px_prev, px);
    // checkLimits(py_prev, py);
    // checkLimits(pz_prev, pz);
    // checkLimits(ry_prev, ry);

    // RCLCPP_INFO(this->get_logger(), "REVISED POSITION----->  X:%f Y %f Z: %f", px, py, pz);
    // RCLCPP_INFO(this->get_logger(), "REVISED GIMBALS----->  X:%f Y %f Z: %f", rx, ry, rz);

    // positionToJointAngles(px, py, pz, dist_sheath_angles_end);

    // // change second option to use gimbals
    // auto rotation = (!dual_sheath_mode) ? 
    // 			(getRotz(dist_sheath_angles_end(1)) * getRoty(dist_sheath_angles_end(0))) : 
    // 			(getRotz(dist_sheath_angles_end(1)) * getRoty(dist_sheath_angles_end(0)));

    // transformToJointVariables(buildTransform(rotation, position), inverse_kinematics_touch_end);

    // RCLCPP_INFO(this->get_logger(), "FEASIBLE: %f TH1D: %f L1: %f TH2D: %f PHI2D: %f L2: %f", 
    // 					inverse_kinematics_touch_end(5),
    // 					inverse_kinematics_touch_end(0), 
    // 					inverse_kinematics_touch_end(1), 
    // 					inverse_kinematics_touch_end(2), 
    // 					inverse_kinematics_touch_end(3), 
    // 					inverse_kinematics_touch_end(4));

    // if (msg.effort[0]) {

    // 	RCLCPP_INFO(this->get_logger(), "GENERATING LINSPACE");
    // 	linspace_step = (inverse_kinematics_touch_end - inverse_kinematics_touch_init) * (1.0 / (LINSPACESTEPS));
    // 	linspace_done = false;
    // }


    //  UNCOMMENT IF DOING BACKUP PLAN
    auto prox_length = 0.0;
    auto prox_theta = 0.0;
    bool run_continuous = true;
    bool dual_sheath_mode = false;
    Eigen::VectorXd dist_sheath_angles_end = Eigen::VectorXd::Zero(6);

    if (run_continuous) {
        // REPLACE NUMBERS HERE WITH HARDCODED NUMBERS FOR PROXIMAL SHEATH
        //if (prox_theta > -1.0 || prox_theta < -48.0) {
        //    prox_theta = prev_prox_ry;
        //}
        //else {
        //    prev_prox_ry = prox_theta;
        //}

        auto rot_use = getRoty(prox_theta);
        Eigen::VectorXd pos_use(4);
        if (dual_sheath_mode) {
            pos_use << (prox_length / prox_theta) * (1.0 - std::cos(deg2rad(prox_theta))), 0, (prox_length / prox_theta)* std::sin(deg2rad(prox_theta)), 1.0;
        }
        else {
            pos_use << 0.0, 0.0, 0.0, 1.0;
        }

        auto transform_use = buildTransform(rot_use, pos_use);
        // RCLCPP_INFO(this->get_logger(), "Received matrix: ");
        // RCLCPP_INFO(this->get_logger(), "%f %f %f %f", 
        // 								transform_use(0,0), transform_use(0,1), transform_use(0,2), transform_use(0,3));
        // RCLCPP_INFO(this->get_logger(), "%f %f %f %f", 
        // 								transform_use(1,0), transform_use(1,1), transform_use(1,2), transform_use(1,3));
        // RCLCPP_INFO(this->get_logger(), "%f %f %f %f", 
        // 								transform_use(2,0), transform_use(2,1), transform_use(2,2), transform_use(2,3));
        // RCLCPP_INFO(this->get_logger(), "%f %f %f %f", 
        // 								transform_use(3,0), transform_use(3,1), transform_use(3,2), transform_use(3,3));

        Eigen::Vector4d cur_pos;
        cur_pos << px, py, pz, 1.0;

        auto new_pos_use = transform_use.inverse() * cur_pos;
        positionToJointAngles(new_pos_use(0), new_pos_use(1), new_pos_use(2), dist_sheath_angles_end);

        dist_sheath_angles_end(4) = prox_length;
        if (dual_sheath_mode) {
            dist_sheath_angles_end(3) = prox_theta;
        }
        else {
            //dist_sheath_angles_end(3) = prev_prox_ry;
        }
    }
    //else {
    //    RCLCPP_INFO(this->get_logger(), "MOVING PROXIMAL SHEATH");
    //    // REPLACE NUMBERS HERE WITH HARDCODED NUMBERS FOR PROXIMAL SHEATH
    //    linspace_step_prox = (20.0 - 0.0) * (1.0 / (LINSPACESTEPS));
    //    linspace_done = false;
    //} */

    /* Eigen::VectorXd q(9);
    double pi = M_PI;
    q << pi / 2, 0, 50, -pi, pi, 20, pi / 4, -3 * pi / 4, 10;
    Eigen::Vector3d pos;
    Eigen::Matrix3d rot;

    forwardKinematics(q, pos, rot);

    //std::cout << pos << std::endl;
    //std::cout << rot << std::endl;

    //std::cout << generateJacobian(q) << std::endl;
    Eigen::VectorXd jointHighLimits(9);
    Eigen::VectorXd jointLowLimits(9);
    Eigen::VectorXd jointMidPoints(9);
    jointHighLimits << q(3-1) / 0.032, 1000 * pi, 0.05, q(6 - 1) / 0.032, 1000 * pi, 0.05, q(9 - 1) / 0.032, 1000 * pi, 0.05;
    jointLowLimits << 0, -1000 * pi, 0.001, -q(6 - 1) / 0.032, -1000 * pi, 0.001, -q(9 - 1) / 0.032, -1000 * pi, 0.001;
    jointMidPoints << pi / 4, 0, 0.025, 0, 0, 0.025, 0, 0, 0.025;

    //std::cout << generateNullSpaceObjectiveVector(q, jointMidPoints, jointHighLimits, jointLowLimits) << std::endl;

    Eigen::Matrix3d g;
    g << 1, 2, 3, 4, 8, 12, 7, 8, 9;
    //std::cout << pinv(g,0.002) << std::endl;

    //std::cout << jointLowLimits.maxCoeff() << std::endl;

    

    Eigen::MatrixXd jacobian = generateJacobian(q);

    //hie(jointMidPoints);

    
    //Eigen::VectorXd ksidesired(6);
    //ksidesired << 0.01, 0.01, 0.01, 0, 0, 0;
    ////std::cout << jacobian << std::endl;
    //Eigen::VectorXd nullspace = generateNullSpaceObjectiveVector(q, jointMidPoints, jointHighLimits, jointLowLimits);
    ////std::cout << nullspace << std::endl;
    //Eigen::VectorXd maxjointspeed(9);
    //maxjointspeed << pi / 6, 1000 * pi, 0.003, pi / 6, 1000 * pi, 0.003, pi / 6, 1000 * pi, 0.003;

    //Eigen::VectorXd velocities = calculateJointVelocities(jacobian, ksidesired, nullspace, maxjointspeed);
    //std::cout << velocities << std::endl;
    //std::cout << "hihihihi" << std::endl;
    //std::cout << velocities.segment(3, 5) << std::endl;

    //double l[3] = {1, 2, 3};
    //double h[3] = { 3, 4, 5 };

    //auto b = h - l;
    //Eigen::ArrayXd l(3);
    //l << 1, 2, 3;
    //Eigen::ArrayXd h(3);
    //h << 3, 4, 5;

    //Eigen::ArrayXd m(3);
    //m << 3, 3, 2;

    ////Eigen::ArrayXd g(3);
    //
    //auto g = (m > (h - l));

    //std::cout << g << std::endl;

    //Eigen::ArrayXd test = Eigen::ArrayXd::Map(jointHighLimits.data(), jointHighLimits.size());
    //std::cout << test << std::endl;

    int numJoints = numTubes * numConfigParams;
    double minAllowedRadii = 0.032;
    Eigen::VectorXd globalJointLowLimits(numJoints); // minimum allowable lengths/angles of the joints
    globalJointLowLimits << 0.0000000000000001, -1000.0 * M_PI, 0.01, -q(5) / minAllowedRadii, -1000.0 * M_PI, 0.01, 0.001, 0.0, 0.001;

    Eigen::VectorXd globalJointHighLimits(numJoints); // maximum allowable lengths/angles of the joints
    globalJointHighLimits << q(2) / minAllowedRadii, 1000.0 * M_PI, 0.05, q(5) / minAllowedRadii, 1000.0 * M_PI, 0.05, 0.001, 0.0, 0.001;
    
    //Eigen::VectorXd exceededJoints = Eigen::VectorXd::Zero(numJoints); // keeps running count of joints that have been checked
    Eigen::VectorXd limitHighJoint = globalJointHighLimits;
    Eigen::VectorXd limitLowJoint = globalJointLowLimits;

    Eigen::VectorXd qNew(9);
    qNew << pi / 2, 0, 50, -pi, pi, 20, pi / 4, -3 * pi / 4, 10;

    for (int i = 0; i < numTubes; ++i) {
        if (i > 0) {
            limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numTubes), -1.0 * (qNew(i * numTubes + 2) / minAllowedRadii)); // theta
        }
        else {
            limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numTubes), 0.0); // theta
        }
        limitLowJoint(i * numConfigParams + 2) = std::max(globalJointLowLimits(i * numTubes + 2), (qNew(i * numTubes) * minAllowedRadii)); // length

        limitHighJoint(i * numConfigParams) = std::min(globalJointHighLimits(i * numTubes), (qNew(i * numTubes + 2) / minAllowedRadii)); // theta
        limitHighJoint(i * numConfigParams + 2) = globalJointHighLimits(i * numTubes + 2); // length
    }

    int exceededJoint = numJoints;
    for (int j = 0; j < numJoints; ++j) {
        if (((qNew(j) - limitLowJoint(j)) / (limitHighJoint(j) - limitLowJoint(j))) > 1) {
            exceededJoint = j;
            break;
        }
    }

    std::cout << exceededJoint << std::endl;

    Eigen::VectorXd a(3);
    a << 2, 2, 2;
    Eigen::VectorXd b(3);
    b << 4, 5, 6;

    Eigen::VectorXd c = Eigen::VectorXd::Ones(3);
    Eigen::VectorXd d(3);
    d = c.array() + a.array() * b.array();

    
    std::cout << d << std::endl; */

    /* int numJoints = numConfigParams * numTubes;
    Eigen::VectorXd q(9);
    q << 1.5625, 1.5, 0.05, M_PI / 4, 1.5, 0.045, 0.001, 0.0, 0.001;

    Eigen::VectorXd globalJointLowLimits(numJoints); // minimum allowable lengths/angles of the joints
    globalJointLowLimits << 0.0000000000000001, -1000.0 * M_PI, 0.001, -q(5) / 0.032, -1000.0 * M_PI, 0.001, 0.001, 0.0, 0.001;

    // these must be hardcoded; CHANGE THESE WHEN CHANGING NUMBER OF ACTIVE TUBES
    Eigen::VectorXd globalJointHighLimits(numJoints); // maximum allowable lengths/angles of the joints
    globalJointHighLimits << q(2) / 0.032, 1000.0 * M_PI, 0.05, q(5) / 0.032, 1000.0 * M_PI, 0.05, 0.001, 0.0, 0.001;

    // perform differential kinematics calculation to obtain Jacobian
    Eigen::MatrixXd J_normal = generateJacobian(q);

    Eigen::MatrixXd J_use;          // Jacobian that will be used for robot control
    Eigen::VectorXd tipVelocityDesired_use(3); // desired input that will be used for robot control

    J_use = J_normal.block(0, 0, 3, numConfigParams * numTubes);

    std::cout << "jacobian: " << J_use << std::endl;
    tipVelocityDesired_use << -0.001, -0.001, -0.001;

    Eigen::VectorXd jointMidPoints(9);
    jointMidPoints << M_PI / 4.0, 0, 0.025, 0.0, 0.0, 0.025, 0.0, 0.0, 0.025;
    // calculate secondary task for null space manipulation (redundancy resolution)
    Eigen::VectorXd nullspaceObjectiveVector = generateNullSpaceObjectiveVector(q,
        jointMidPoints,
        globalJointHighLimits,
        globalJointLowLimits);

    std::cout << "redundancy resolution: " << nullspaceObjectiveVector << std::endl;

     // this must be hardcoded; CHANGE THIS WHEN CHANGING NUMBER OF ACTIVE TUBES
    for (int ii : {1, 6, 7, 8}) {
        J_use.col(ii) = Eigen::VectorXd::Zero(J_use.rows());
        nullspaceObjectiveVector(ii) = 0.0;
    }

    // calculate joint velocities
    Eigen::VectorXd qDot = calculateJointVelocities(J_use,
        tipVelocityDesired_use,
        nullspaceObjectiveVector);

    Eigen::VectorXd qNew = q + qDot * 0.01; // recalculate config params
    bool applyJointLimits = true;
    int counter = 0;

    std::cout << "qnew is: " << qNew << std::endl;

    if (applyJointLimits) {
        bool finishedChecking = false;
        //int numJoints = numTubes * numConfigParams;
        Eigen::VectorXd limitHighJoint = globalJointHighLimits;
        Eigen::VectorXd limitLowJoint = globalJointLowLimits;

        while (!finishedChecking) {

            // calculate joint limits dynamically
            for (int i = 0; i < numTubes; ++i) {
                if (i > 0) {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numTubes), -1.0 * (q(i * numTubes + 2) / 0.032)); // theta
                }
                else {
                    limitLowJoint(i * numConfigParams) = std::max(globalJointLowLimits(i * numTubes), 0.0); // theta
                }
                limitLowJoint(i * numConfigParams + 2) = std::max(globalJointLowLimits(i * numTubes + 2), (q(i * numTubes) * 0.032)); // length

                limitHighJoint(i * numConfigParams) = std::min(globalJointHighLimits(i * numTubes), (q(i * numTubes + 2) / 0.032)); // theta
                limitHighJoint(i * numConfigParams + 2) = globalJointHighLimits(i * numTubes + 2); // length
            }

            std::cout << "low limits are " << limitLowJoint << std::endl;

            std::cout << "high limits joints are: " << limitHighJoint << std::endl;

            int exceededJoint = numJoints;
            // find and save first joint that exceeds previously calculated joint limits
            for (int j = 0; j < numJoints; ++j) {
                if ((std::abs(qNew(j) - limitLowJoint(j)) / std::abs(limitHighJoint(j) - limitLowJoint(j))) > 1) {
                    exceededJoint = j;
                    std::cout << "exceeded joint in for loop is: " << (j + 1) << std::endl;
                    break;
                }
            }

            if (exceededJoint != numJoints) { // if a joint has exceeded the limits
                J_use.col(exceededJoint) = Eigen::VectorXd::Zero(J_use.rows()); // wipe the column in the Jacobian of the current param
                nullspaceObjectiveVector(exceededJoint) = 0.0; // wipe q_0_dot

                std::cout << "joint exceeding limits: " << J_use << std::endl;

                // calculate joint velocities
                qDot = calculateJointVelocities(J_use,
                    tipVelocityDesired_use,
                    nullspaceObjectiveVector);

                qNew = q + qDot * 0.01; // recalculate config params

                // std::string minOrMax = (qNewPre(exceededJoint) < limitLowJoint(exceededJoint)) ? "MIN" : "MAX";
                // RCLCPP_INFO(this->get_logger(), "THETA %s LIMIT REACHED FOR SEGMENT %i", minOrMax.c_str(), exceededJoint / 3);
                // RCLCPP_INFO(this->get_logger(), "THETA VALUE %i IS: %f", exceededJoint, qNewPre(exceededJoint));
            }
            else {
                finishedChecking = true;
            }

            counter++;
        }
        std::cout << "num times going through loop is: " << counter << std::endl;
    } */

    /* Eigen::VectorXd q(5);
    q << 1, 2, 3, 4, 5;
    std::cout << q.segment(0, 2) << std::endl;

    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> capitalPiMatrixi_[3];
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> invCapitalPiMatrixi_[3];

    Eigen::MatrixXd beta(numTotalTubes, numTendonsPerTube);
    // this must be hardcoded
    beta.row(0) << M_PI/6.0, M_PI/4.0;                // outer sheath
    beta.row(1) << M_PI / 4.0, 3.0 * M_PI / 4.0;          // intermediate sheath
    beta.row(2) << 0.0, M_PI / 2.0;                       // TODO: inner sheath (currently unused)

    // this must be hardcoded
    for (int i = 0; i < numTotalTubes; ++i)
    {
        capitalPiMatrixi_[i] << std::cos(beta(i, 0)),
            std::cos(beta(i, 1)),
            std::sin(beta(i, 0)),
            std::sin(beta(i, 1));

        invCapitalPiMatrixi_[i] = capitalPiMatrixi_[i].inverse();
    }

    Eigen::VectorXd w(2);
    w << 1, 2;

    std::cout << invCapitalPiMatrixi_[0] << std::endl;

    std::cout << invCapitalPiMatrixi_[0] * w << std::endl; */

    /* //Eigen::VectorXd a(10);
    //a << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    //auto b = 2 * a;
    //
    //time_t tstart, tend;
    //tstart = time(0);

    //int j = 730;
    //for (int i = 0; i < 1000; ++i) {
    //    if (i == j) {
    //        tend = time(0);
    //        std::cout << "hi and time taken is: " << difftime(tend, tstart) << std::endl;
    //    }
    //}

    Eigen::VectorXd tendonTotalDisplacement_(12);
    tendonTotalDisplacement_ << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    Eigen::VectorXd length(3);
    length << 1, 2, 3;
    Eigen::VectorXd a(3);
    Eigen::VectorXd b(3);

    //for (int i = 0; i < 8; ++i) {
    //    auto tosend = tendonTotalDisplacement_ * (i * 0.005);
    //    auto tosend2 = length * (i * 0.01);
    //    jointSpaceToConfigurationSpace(3, tosend, tosend2, a, b);
    //}
    for (int k = 0; k < numTubes; ++k) {
        std::cout << tendonTotalDisplacement_.segment(k * numTendonsPerTube + 3, numTendonsPerTube) << std::endl;
    }

    std::vector<int> g;
    g.push_back(2);
    g.push_back(5);
    g.push_back(6);
    g.push_back(8);

    for (int i : g) {
        std::cout << tendonTotalDisplacement_(i) << std::endl;
    }

    std::cout << "rem: " << 5 % 3 << std::endl;

    if (std::find(g.begin(), g.end(), 2) != g.end()) {
        std::cout << "hi" << std::endl;
    }
    else {
        std::cout << "hello" << std::endl;
    } */

    fesetround(FE_TONEAREST); // Round to nearest (even) like MATLAB

    theta_init_ << deg2rad(25), deg2rad(25), 0.00005; // initial bending angles of each sheath
    phi_init_ << 180.0 * M_PI / 180.0, -90.0 * M_PI / 180.0, 0.0 * M_PI / 180.0; // initial bending plane offsets of each sheath
    length_init_ << 10.0, 10.0, 10.0; // initial length of the steerable section of each sheath in mm
    delta_ << 2.3, 1.8, 1.0; // radius (in meters) of cross sectional circle centered at sheath centerline and passes through the tendon locations

    Eigen::VectorXd q(numConfigParams * numTubes);
    Eigen::VectorXd qDot(numConfigParams * numTubes);

    jointMidPoints << M_PI / 4.0, 0, 25.0, 0.0, 0.0, 25.0, 0.0, 0.0, 25.0;

    maxJointSpeed << M_PI / 6.0, 1000.0 * M_PI, 3,
        M_PI / 6.0, 1000.0 * M_PI, 3,
        M_PI / 6.0, 1000.0 * M_PI, 3;

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

    for (int i = 0; i < numTubes; ++i) {
        q_sim.segment((i * numConfigParams), numConfigParams) << theta_init_(i), phi_init_(i), length_init_(i);
    }

    std::vector<MedianFilter> medfilt;
    for (int i = 0; i < 6; ++i) {
        medfilt.push_back(MedianFilter(100));
    }

    std::cout << q_sim << std::endl;
    Eigen::Matrix3d dummyRotMM_i_wrtG1 = Eigen::Matrix3d::Identity();

    //Eigen::VectorXd q_test(9);
    //q_test << deg2rad(67), deg2rad(54), 20, deg2rad(67), deg2rad(54), 20, deg2rad(67), deg2rad(54), 20;

    forwardKinematics(q_sim, pTip_init_sim, dummyRotMM_i_wrtG1);
    //std::cout << "ptip: " << pTip_init_sim << std::endl;
    //std::cout << "rot: " << dummyRotMM_i_wrtG1 << std::endl;

    //std::cout << generateJacobian(q_sim) << std::endl;

    std::vector<Eigen::VectorXd> printqsims;
    std::vector<Eigen::VectorXd> printqdotsims;

    for (int j = 0; j < 10000; ++j) {
        auto time = j * deltaT;
        Eigen::Matrix3d dummyRotMM_i_wrtG = Eigen::Matrix3d::Identity();
        
        forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);

        std::pair<Eigen::Vector3d, Eigen::Vector3d> result =
            calculate_circleTrajectory_position_velocity_f(pTip_init_sim, time, 0.25, outOfWorkspace, pTip_sim);
        automated_velocity = result.second;

        Eigen::VectorXd tipVelocityDesired(6); // command from manual control
        tipVelocityDesired << automated_velocity(0), automated_velocity(1), automated_velocity(2), 0.0, 0.35 * (2.0 * M_PI / 500.0) * std::cos(j * ((2.0 * M_PI) / 500.0)), 0.0;

        Eigen::VectorXd displacementCalc(nOfDrivenMotors);
        displacementCalc = Eigen::VectorXd::Zero(nOfDrivenMotors); // motor displacement

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
        Eigen::MatrixXd jacobian = generateJacobian(q_sim);


        forwardKinematics(q_sim, pTip_sim, dummyRotMM_i_wrtG);
        writeVectorToCSV("C:/Users/ch256744/BCH Dropbox/Phillip Tran/ExtendedICRAPaper/AbdulsCode/MatlabCode/endtip_precise1021_20circle_2.csv", pTip_sim);
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

    Mat m(1, 1);
    long double u = 5;
    m(0, 0) = u;
    std::cout << m << std::endl;

    //std::vector<double> locked_joints;
    //locked_joints.push_back(1);
    //locked_joints.push_back(6);
    //locked_joints.push_back(7);
    //locked_joints.push_back(8);
    //

    //std::cout << find_element(locked_joints, 1) << std::endl; */

    /* //double f = 94;  // Known value of f
    //double a = 1.0;  // Known value of a
    //double b = 2.0;  // Known value of b
    //double c = 3.0;  // Known value of c

    Timer timer;

    //double th = solve_equation(f, a, b, c);

    //if (!std::isnan(th)) {
    //    std::cout << "Solution: th = " << th << " radians" << std::endl;
    //    std::cout << "th = " << th * 180 / PI << " degrees" << std::endl;

    //    // Verify the solution
    //    double result = a + b * std::atan(th) + c * th;
    //    std::cout << "Verification: f = " << result << std::endl;
    //}

    double f1_target = 47.04, f2_target = 109.65;
    double a1 = 1.0, b1 = 2.0, c1 = 3.0;
    double a2 = 2.0, b2 = 3.0, c2 = 4.0;

    timer.tic();
    std::vector<double> solution = solveEquations(f1_target, f2_target, a1, b1, c1, a2, b2, c2);
    double elapsed_time = timer.toc();

    std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;

    if (!solution.empty()) {
        std::cout << "Solution found:" << std::endl;
        std::cout << "th = " << solution[0] << std::endl;
        std::cout << "phi = " << solution[1] << std::endl;
    }
    else {
        std::cout << "No solution found within the specified iterations." << std::endl;
    } */

    /* 
    //for (int i = 0; i < 1; ++i) {
    //    for (int j = 0; j < 1; ++j) {
    //        Eigen::Matrix3d rot = getRoty(double(j)) * getRotx(double(i)) * getRotz(0.0);
    //        auto c = rotationMatrixToEulerAngles(rot);
    //        std::cout << rad2deg(c(0)) << rad2deg(c(1)) << rad2deg(c(2)) << std::endl;
    //    }
    //}

    Eigen::Matrix3d f;
    f << 0.8763,    0.3031,    0.3745,
        0.3563,    0.1158, - 0.9272,
        - 0.3244,    0.9459, - 0.0065;

    //std::cout << f << std::endl;
    //auto b = rotationToEuler_yxz(f,0,deg2rad(90.8),0);
    //std::cout << rad2deg(b(0)) << " " << rad2deg(b(1)) << " " << rad2deg(b(2)) << std::endl;

    Eigen::MatrixXd y(6,6);
    y << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36;
    auto J_use = y.block(0, 0, 3, 6);
    auto other = y.block(4, 0, 1, 6);

    std::cout << y << std::endl;
    std::cout << other << std::endl;

    Eigen::MatrixXd b(J_use.rows() + other.rows(), J_use.cols());
    b.topRows(J_use.rows()) = J_use;
    b.bottomRows(other.rows()) = other;

    std::cout << b << std::endl;

    std::cout << sgn(0) << std::endl;
    std::cout << sgn(-3) << std::endl; */
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
