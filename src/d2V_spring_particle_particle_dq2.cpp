#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::Matrix3d He;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double squared_dist = (q0-q1).squaredNorm();
    double dist = (q0-q1).norm();
    Eigen::Vector3d q01=q0-q1;
    He = stiffness*q01*q01.transpose()/squared_dist + stiffness*(1.0 - l0/dist)*(I-q01*q01.transpose()/squared_dist);
    H.block(0,0,3,3)=He;
    H.block(0,3,3,3)=-He;
    H.block(3,0,3,3) = -He;
    H.block(3,3,3,3)=He;
    H=-1.0*H;
    
}