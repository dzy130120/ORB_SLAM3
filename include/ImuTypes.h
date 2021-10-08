/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2020 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef IMUTYPES_H
#define IMUTYPES_H

#include<vector>
#include<utility>
#include<opencv2/core/core.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <mutex>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace ORB_SLAM3
{

namespace IMU
{

const float GRAVITY_VALUE=9.81;

//IMU measurement (gyro, accelerometer and timestamp)
//imu数据点
class Point//imu数据点
{
public:
    Point(const float &acc_x, const float &acc_y, const float &acc_z,
             const float &ang_vel_x, const float &ang_vel_y, const float &ang_vel_z,
             const double &timestamp): a(acc_x,acc_y,acc_z), w(ang_vel_x,ang_vel_y,ang_vel_z), t(timestamp){}
    Point(const cv::Point3f Acc, const cv::Point3f Gyro, const double &timestamp):
        a(Acc.x,Acc.y,Acc.z), w(Gyro.x,Gyro.y,Gyro.z), t(timestamp){}
public:
    cv::Point3f a;
    cv::Point3f w;
    double t;
};

//IMU biases (gyro and accelerometer)
//IMU 偏置
class Bias
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & bax;
        ar & bay;
        ar & baz;

        ar & bwx;
        ar & bwy;
        ar & bwz;
    }

public:
    Bias():bax(0),bay(0),baz(0),bwx(0),bwy(0),bwz(0){}
    Bias(const float &b_acc_x, const float &b_acc_y, const float &b_acc_z,
            const float &b_ang_vel_x, const float &b_ang_vel_y, const float &b_ang_vel_z):
            bax(b_acc_x), bay(b_acc_y), baz(b_acc_z), bwx(b_ang_vel_x), bwy(b_ang_vel_y), bwz(b_ang_vel_z){}
    void CopyFrom(Bias &b);
    friend std::ostream& operator<< (std::ostream &out, const Bias &b);

public:
    float bax, bay, baz;
    float bwx, bwy, bwz;
};

//IMU calibration (Tbc, Tcb, noise)
//IMU与camera外参
class Calib
{
    template<class Archive>
    void serializeMatrix(Archive &ar, cv::Mat& mat, const unsigned int version)
    {
        int cols, rows, type;
        bool continuous;

        if (Archive::is_saving::value) {
            cols = mat.cols; rows = mat.rows; type = mat.type();
            continuous = mat.isContinuous();
        }

        ar & cols & rows & type & continuous;
        if (Archive::is_loading::value)
            mat.create(rows, cols, type);

        if (continuous) {
            const unsigned int data_size = rows * cols * mat.elemSize();
            ar & boost::serialization::make_array(mat.ptr(), data_size);
        } else {
            const unsigned int row_size = cols*mat.elemSize();
            for (int i = 0; i < rows; i++) {
                ar & boost::serialization::make_array(mat.ptr(i), row_size);
            }
        }
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        serializeMatrix(ar,Tcb,version);
        serializeMatrix(ar,Tbc,version);
        serializeMatrix(ar,Cov,version);
        serializeMatrix(ar,CovWalk,version);
    }

public:
    Calib(const cv::Mat &Tbc_, const float &ng, const float &na, const float &ngw, const float &naw)
    {
        Set(Tbc_,ng,na,ngw,naw);
    }
    Calib(const Calib &calib);
    Calib(){}

    void Set(const cv::Mat &Tbc_, const float &ng, const float &na, const float &ngw, const float &naw);

public:
    cv::Mat Tcb;
    cv::Mat Tbc;
    cv::Mat Cov, CovWalk;
};

//Integration of 1 gyro measurement
//角速度积分
class IntegratedRotation
{
public:
    IntegratedRotation(){}
    IntegratedRotation(const cv::Point3f &angVel, const Bias &imuBias, const float &time);

public:
    float deltaT; //integration time
    cv::Mat deltaR; //integrated rotation
    cv::Mat rightJ; // right jacobian
};

//Preintegration of Imu Measurements
class Preintegrated
{
    template<class Archive>
    void serializeMatrix(Archive &ar, cv::Mat& mat, const unsigned int version)
    {
        int cols, rows, type;
        bool continuous;

        if (Archive::is_saving::value) {
            cols = mat.cols; rows = mat.rows; type = mat.type();
            continuous = mat.isContinuous();
        }

        ar & cols & rows & type & continuous;
        if (Archive::is_loading::value)
            mat.create(rows, cols, type);

        if (continuous) {
            const unsigned int data_size = rows * cols * mat.elemSize();
            ar & boost::serialization::make_array(mat.ptr(), data_size);
        } else {
            const unsigned int row_size = cols*mat.elemSize();
            for (int i = 0; i < rows; i++) {
                ar & boost::serialization::make_array(mat.ptr(i), row_size);
            }
        }
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & dT;
        serializeMatrix(ar,C,version);
        serializeMatrix(ar,Info,version);
        serializeMatrix(ar,Nga,version);
        serializeMatrix(ar,NgaWalk,version);
        ar & b;
        serializeMatrix(ar,dR,version);
        serializeMatrix(ar,dV,version);
        serializeMatrix(ar,dP,version);
        serializeMatrix(ar,JRg,version);
        serializeMatrix(ar,JVg,version);
        serializeMatrix(ar,JVa,version);
        serializeMatrix(ar,JPg,version);
        serializeMatrix(ar,JPa,version);
        serializeMatrix(ar,avgA,version);
        serializeMatrix(ar,avgW,version);

        ar & bu;
        serializeMatrix(ar,db,version);
        ar & mvMeasurements;
    }

public:
    Preintegrated(const Bias &b_, const Calib &calib);//构造
    Preintegrated(Preintegrated* pImuPre);//拷贝构造
    Preintegrated() {}//默认构造
    ~Preintegrated() {}//析构
    void CopyFrom(Preintegrated* pImuPre);
    void Initialize(const Bias &b_);
    void IntegrateNewMeasurement(const cv::Point3f &acceleration, const cv::Point3f &angVel, const float &dt);
    void Reintegrate();
    void MergePrevious(Preintegrated* pPrev);
    void SetNewBias(const Bias &bu_);//设置新的bias
    IMU::Bias GetDeltaBias(const Bias &b_);
    cv::Mat GetDeltaRotation(const Bias &b_);
    cv::Mat GetDeltaVelocity(const Bias &b_);
    cv::Mat GetDeltaPosition(const Bias &b_);
    cv::Mat GetUpdatedDeltaRotation();
    cv::Mat GetUpdatedDeltaVelocity();
    cv::Mat GetUpdatedDeltaPosition();
    cv::Mat GetOriginalDeltaRotation();
    cv::Mat GetOriginalDeltaVelocity();
    cv::Mat GetOriginalDeltaPosition();
    Eigen::Matrix<double,15,15> GetInformationMatrix();
    cv::Mat GetDeltaBias();
    Bias GetOriginalBias();
    Bias GetUpdatedBias();

public:
    float dT;
    cv::Mat C; //协方差矩阵
    cv::Mat Info;//信息矩阵
    cv::Mat Nga, NgaWalk;//噪声、随机游走

    // Values for the original bias (when integration was computed)
    Bias b;//更新前bias
    cv::Mat dR, dV, dP;//预积分
    cv::Mat JRg, JVg, JVa, JPg, JPa;
    cv::Mat avgA;//平均加速度
    cv::Mat avgW;//平均角速度


private:
    // Updated bias
    Bias bu;//更新后bias
    // Dif between original and updated bias
    // This is used to compute the updated values of the preintegration
    cv::Mat db;//bias更新前后的变化量

    struct integrable
    {
        integrable(const cv::Point3f &a_, const cv::Point3f &w_ , const float &t_):a(a_),w(w_),t(t_){}
        cv::Point3f a;
        cv::Point3f w;
        float t;
    };

    std::vector<integrable> mvMeasurements;//imu多组读数

    std::mutex mMutex;
};

// Lie Algebra Functions
cv::Mat ExpSO3(const float &x, const float &y, const float &z);
Eigen::Matrix<double,3,3> ExpSO3(const double &x, const double &y, const double &z);
cv::Mat ExpSO3(const cv::Mat &v);
cv::Mat LogSO3(const cv::Mat &R);
cv::Mat RightJacobianSO3(const float &x, const float &y, const float &z);
cv::Mat RightJacobianSO3(const cv::Mat &v);
cv::Mat InverseRightJacobianSO3(const float &x, const float &y, const float &z);
cv::Mat InverseRightJacobianSO3(const cv::Mat &v);
cv::Mat Skew(const cv::Mat &v);
cv::Mat NormalizeRotation(const cv::Mat &R);

}

} //namespace ORB_SLAM2

#endif // IMUTYPES_H
