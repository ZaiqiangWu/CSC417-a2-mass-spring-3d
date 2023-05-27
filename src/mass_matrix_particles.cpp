#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    std::vector<Eigen::Triplet<double>> tripletlist;
    int n = q.rows();
    M.resize(n,n);
    for(int i=0;i<n;++i)
    {
        tripletlist.emplace_back(i,i,mass);
    }
    M.setFromTriplets(tripletlist.begin(),tripletlist.end());
    M.makeCompressed();

    
}
