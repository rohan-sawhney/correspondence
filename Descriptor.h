#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include "Types.h"
#define HKS 0
#define FAST_HKS 1
#define WKS 2

class Descriptor {
public:
    // constructor
    Descriptor();
    
    // setup
    void setup(Mesh *mesh0);
    
    // compute
    void compute(int descriptor);
    
private:
    // compute eigenvalues and eigenvectors
    void computeEig(const Eigen::SparseMatrix<double>& W,
                    const Eigen::SparseMatrix<double>& A);
    
    // compute hks
    void computeHks();
    
    // compute fast hks
    void computeFastHks();
    
    // compute wks
    void computeWks();
    
    // normalize
    void normalize();
    
    // member variable
    Mesh *mesh;
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
};

#endif
