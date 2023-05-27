#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        int n = l0.rows();
        f.resize(q.rows());
        f.setZero();
        for(int e=0;e<n;++e)
        {
            Eigen::Vector6d local_f;
            int id0=E(e,0);
            int id1=E(e,1);
            dV_spring_particle_particle_dq(local_f,q.block(3*id0,0,3,1),q.block(3*id1,0,3,1),l0(e),k);
            f.block(3*id0,0,3,1)+=local_f.block(0,0,3,1);
            f.block(3*id1,0,3,1)+=local_f.block(3,0,3,1);
        }
        int n_vert = q.rows()/3;
        for(int v=0;v<n_vert;++v)
        {
            Eigen::Vector3d gravity;
            dV_gravity_particle_dq(gravity,mass,Eigen::Vector3d(0,0,-10.0));
            f.block(3*v,0,3,1)+=gravity;
        }



    };