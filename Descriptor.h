#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include "Types.h"
#define HKS 0
#define WKS 1

class Descriptor {
public:
    // constructor
    Descriptor();
    
    // setup
    void setup(Mesh *mesh0, int k0);
    
    // compute
    void compute(int n, int descriptorName);
    
private:
    // build laplace operator
    void buildLaplace(Eigen::SparseMatrix<double>& L);
    
    // build area matrix
    void buildAreaMatrix(Eigen::SparseMatrix<double>& A);
    
    // compute eigenvalues and eigenvectors
    void computeEig(const Eigen::SparseMatrix<double>& L,
                    const Eigen::SparseMatrix<double>& A);
    
    // compute hks
    void computeHks(int n);
    
    // compute wks
    void computeWks(int n);
    
    // normalize
    void normalize();
    
    // member variable
    Mesh *mesh;
    int k;
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
};

#endif
