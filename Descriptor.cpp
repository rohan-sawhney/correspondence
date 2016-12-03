#include "Descriptor.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "MultiresMesh.h"
#include <spectra/include/MatOp/SparseSymMatProd.h>
#include <spectra/include/MatOp/SparseCholesky.h>
#include <spectra/include/SymGEigsSolver.h>
#define K 100
#define N 10

Descriptor::Descriptor()
{
    
}

void buildLaplace(Mesh *mesh, Eigen::SparseMatrix<double>& L)
{
    std::vector<Eigen::Triplet<double>> LTriplets;
    
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        do {
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan());
            sumCoefficients += coefficient;
            
            LTriplets.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -coefficient));
            
            he = he->flip->next;
        } while (he != v->he);
        
        LTriplets.push_back(Eigen::Triplet<double>(v->index, v->index, sumCoefficients));
    }
    
    L.setFromTriplets(LTriplets.begin(), LTriplets.end());
}

void buildAreaMatrix(Mesh *mesh, Eigen::SparseMatrix<double>& A)
{
    std::vector<Eigen::Triplet<double>> ATriplets;
    
    double sum = 0.0;
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        double area = v->dualArea();
        ATriplets.push_back(Eigen::Triplet<double>(v->index, v->index, area));
        sum += area;
    }
    
    A.setFromTriplets(ATriplets.begin(), ATriplets.end());
    A *= mesh->vertices.size() / sum;
}

void Descriptor::computeEig(const Eigen::SparseMatrix<double>& L,
                            const Eigen::SparseMatrix<double>& A)
{
    Spectra::SparseSymMatProd<double> opL(L);
    Spectra::SparseCholesky<double> opA(A);
    
    Spectra::SymGEigsSolver<double,
    Spectra::SMALLEST_MAGN,
    Spectra::SparseSymMatProd<double>,
    Spectra::SparseCholesky<double>,
    Spectra::GEIGS_CHOLESKY> geigs(&opL, &opA, K, 2*K);
    
    geigs.init();
    geigs.compute();
    
    if (geigs.info() == Spectra::SUCCESSFUL) {
        evals = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }
    
    Eigen::MatrixXd err = L*evecs - A*evecs*evals.asDiagonal();
    std::cout << "||Lx - λAx||_inf = " << err.array().abs().maxCoeff() << std::endl;
}

void Descriptor::setup(Mesh *mesh0)
{
    mesh = mesh0;
    
    std::string filename = mesh->name;
    filename.replace(filename.find_last_of(".")+1, 3, "eig");
    std::ifstream in(filename);
    
    if (in.is_open()) {
        // read eigvalues and eigenvectors from file
        MeshIO::readEig(in, evals, evecs);
        
    } else {
        int v = (int)mesh->vertices.size();
        
        // build laplace operator
        Eigen::SparseMatrix<double> L(v, v);
        buildLaplace(mesh, L);
        
        // build area matrix
        Eigen::SparseMatrix<double> A(v, v);
        buildAreaMatrix(mesh, A);
        
        // compute eigenvectors and eigenvalues
        computeEig(L, A);
        
        // write eigvalues and eigenvectors to file
        std::ofstream out(filename);
        if (out.is_open()) MeshIO::writeEig(out, evals, evecs);
        out.close();
    }
    
    in.close();
}

void Descriptor::computeHks()
{
    const double ln = 4*log(10);
    const double tmin = ln/evals(0);
    const double step = (ln/evals(K-2) - tmin) / N;
    
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->descriptor = Eigen::VectorXd::Zero(N);
        Eigen::VectorXd C = Eigen::VectorXd::Zero(N);
        
        for (int j = K-2; j >= 0; j--) {
            double phi2 = evecs(v->index, j)*evecs(v->index, j);
            double t = tmin;
            double factor = 0.5;
            
            for (int i = 0; i < N; i++) {
                double exponent = exp(-evals(j)*t);
                v->descriptor(i) += phi2*exponent;
                C(i) += exponent;
                t += factor*step;
                
                // take larger steps with increasing t to bias ts towards high frequency features
                factor += 0.1;
            }
        }
        
        // normalize
        for (int i = 0; i < N; i++) {
            v->descriptor(i) /= C(i);
        }
    }
}

void Descriptor::extrapolateEvals(double& xhat, double& yhat, double& m)
{
    // compute averages
    xhat = 0.0; yhat = 0.0;
    for (int i = 0; i < K; i++) {
        xhat += i;
        yhat += evals(i);
    }
    xhat /= K; yhat /= K;
    
    // compute slope
    double den = 0.0; m = 0.0;
    for (int i = 0; i < K; i++) {
        m += (i - xhat)*(evals(i) - yhat);
        den += (i - xhat)*(i - xhat);
    }
    m /= den;
}

void Descriptor::computeFastHks()
{
    const int d = 2;
    const int C = 1000;
    const double c = 0.2;
    const double reps = 1e-4;
    const std::vector<int> t = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    
    // extrapolate eigenvalues
    double xhat, yhat, m;
    extrapolateEvals(xhat, yhat, m);
    
    // build mutliresolution structure
    MultiresMesh mrm(mesh, d, C);
    mrm.build();
    
    for (int i = 0; i < (int)t.size(); i++) {
        // choose coarsest resolution level l s.t. cn > r(t)
        int r = 0, l = mrm.numLods()-1;
        while (exp(-(yhat - m*(r - xhat))*t[i]) > reps) r++;
        while (l > 0 && c*mrm.lod(l)->vertices.size() < r) l--;
        
        // TODO: compute sparse heat kernel on h
        
        // TODO: project sparse heat kernel on finest resolution level
    }
}

void Descriptor::computeWks()
{
    Eigen::VectorXd logE(K);
    for (int i = 0; i < K; i++) logE(i) = log(evals(i));
    const double emin = logE(K-2);
    const double emax = logE(0)/1.02;
    const double step = (emax - emin) / N;
    const double sigma22 = 2*pow(7*(emax - emin) / K, 2);
 
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->descriptor = Eigen::VectorXd::Zero(N);
        Eigen::VectorXd C = Eigen::VectorXd::Zero(N);

        for (int j = 0; j < K-1; j++) {
            double phi2 = evecs(v->index, j)*evecs(v->index, j);
            double e = emin;
            double factor = 1.5;
            
            for (int i = 0; i < N; i++) {
                double exponent = exp(-pow(e - logE(j), 2) / sigma22);
                v->descriptor(i) += phi2*exponent;
                C(i) += exponent;
                e += factor*step;
                
                // take smaller steps with increasing e to bias es towards high frequency features
                factor -= 0.1;
            }
        }
        
        // normalize
        for (int i = 0; i < N; i++) {
            v->descriptor(i) /= C(i);
        }
    }
}

void Descriptor::normalize()
{
    int n = (int)mesh->vertices[0].descriptor.size();
    for (int i = 0; i < n; i++) {
        // compute min and max
        double min = 0.0, max = 0.0;
        for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
            min = std::min(min, v->descriptor(i));
            max = std::max(max, v->descriptor(i));
        }
        
        // normalize
        for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
            v->descriptor(i) = (v->descriptor(i) - min) / (max - min);
        }
    }
}

void Descriptor::compute(int descriptorName)
{
    // compute descriptor
    if (descriptorName == HKS) computeHks();
    else if (descriptorName == FAST_HKS) computeFastHks();
    else if (descriptorName == WKS) computeWks();
    
    // normalize
    normalize();
}
