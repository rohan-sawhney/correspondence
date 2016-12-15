#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include "Types.h"
#define HKS 0
#define FAST_HKS 1
#define WKS 2
#define CURVE 3

class Descriptor {
public:
    // constructor
    Descriptor(Mesh *mesh0);
    
    // compute
	void compute(int descriptor, bool loadEig = true, std::string outFilename = "");
    
private:
    // compute eigenvalues and eigenvectors
    void computeEig(int K, bool loadEig);
    
    // compute hks
    void computeHks();
    
    // compute fast hks
    void computeFastHks();
    
    // compute wks
    void computeWks();
    
    // compute curvature descriptor
    void computeCurve();
    
    // normalize
    void normalize();
    
    // member variable
    Mesh *mesh;
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
};

#endif
