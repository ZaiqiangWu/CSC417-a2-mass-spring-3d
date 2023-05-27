#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
        int n_vert = q.rows()/3;
        K.resize(3*n_vert,3*n_vert);

    std::vector<Eigen::Triplet<double>> tripletlist;

    int n_edge = E.rows();
    for (int e=0;e<n_edge;++e)
    {
        Eigen::Matrix66d H;
        int id0=E(e,0);
        int id1 = E(e,1);
        int ids[2]={id0,id1};

        d2V_spring_particle_particle_dq2(H,q.block(3*id0,0,3,1),q.block(3*id1,0,3,1),l0(e),k);
        for(int i=0;i<2;++i)
        {
            for(int j=0;j<2;++j)
            {
                for(int r=0;r<3;++r)
                {
                    for(int c=0;c<3;++c)
                    {
                        tripletlist.emplace_back(3*ids[i]+r,3*ids[j]+c,H(3*i+r,3*j+c));
                    }
                }
            }
        }


    }

    K.setFromTriplets(tripletlist.begin(),tripletlist.end());
    K.makeCompressed();

        
    };