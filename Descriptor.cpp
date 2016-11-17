#include "Descriptor.h"
#include "Mesh.h"
#include "MeshIO.h"
#include <spectra/include/MatOp/SparseSymMatProd.h>
#include <spectra/include/MatOp/SparseCholesky.h>
#include <spectra/include/SymGEigsSolver.h>

Descriptor::Descriptor()
{
    
}

void Descriptor::buildLaplace(Eigen::SparseMatrix<double>& L)
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

void Descriptor::buildAreaMatrix(Eigen::SparseMatrix<double>& A)
{
    std::vector<Eigen::Triplet<double>> ATriplets;
    
    double sum = 0.0;
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        double area = v->dualArea();
        ATriplets.push_back(Eigen::Triplet<double>(v->index, v->index, area));
        sum += area;
    }
    
    A.setFromTriplets(ATriplets.begin(), ATriplets.end());
    A /= sum;
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
    Spectra::GEIGS_CHOLESKY> geigs(&opL, &opA, k, 2*k);
    
    geigs.init();
    geigs.compute();
    
    if (geigs.info() == Spectra::SUCCESSFUL) {
        evals = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }
    
    Eigen::MatrixXd err = L*evecs - A*evecs*evals.asDiagonal();
    std::cout << "||Lx - λAx||_inf = " << err.array().abs().maxCoeff() << std::endl;
}

void Descriptor::setup(Mesh *mesh0, int k0)
{
    mesh = mesh0;
    k = k0;
    
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
        buildLaplace(L);
        
        // build area matrix
        Eigen::SparseMatrix<double> A(v, v);
        buildAreaMatrix(A);
        
        // compute eigenvectors and eigenvalues
        computeEig(L, A);
        
        // write eigvalues and eigenvectors to file
        std::ofstream out(filename);
        if (out.is_open()) MeshIO::writeEig(out, evals, evecs);
        out.close();
    }
    
    in.close();
}

void Descriptor::computeHks(int n)
{
    double ln = 4*log(10);
    double tmin = ln/evals(0);
    double step = (ln/evals(k-2) - tmin) / n;
    
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->feature = Eigen::VectorXd::Zero(n);
        
        for (int j = k-2; j >= 0; j--) {
            double phi2 = evecs(v->index, j)*evecs(v->index, j);
            double t = tmin;
            double factor = 0.5;
            
            for (int i = 0; i < n; i++) {
                v->feature(i) += phi2*exp(-evals(j)*t);
                t += factor*step;
                
                // take larger steps with increasing t to bias ts towards high frequency features
                factor += 0.1;
            }
        }
    }
}

void Descriptor::computeWks(int n)
{
    Eigen::VectorXd logE(k);
    for (int i = 0; i < k; i++) logE(i) = log(evals(i));
    double emin = logE(k-2);
    double emax = logE(0)/1.02;
    double step = (emax - emin) / n;
    double sigma22 = 2*pow(7*(emax - emin) / k, 2);
 
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->feature = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd C = Eigen::VectorXd::Zero(n);

        for (int j = 0; j < k-1; j++) {
            double phi2 = evecs(v->index, j)*evecs(v->index, j);
            double e = emin;
            double factor = 1.5;
            
            for (int i = 0; i < n; i++) {
                double exponent = exp(-pow(e - logE(j), 2) / sigma22);
                v->feature(i) += phi2*exponent;
                C(i) += exponent;
                e += factor*step;
                
                // take smaller steps with increasing e to bias es towards high frequency features
                factor -= 0.1;
            }
        }
        
        // normalize
        for (int i = 0; i < n; i++) {
            v->feature(i) /= C(i);
        }
    }
}

void Descriptor::normalize()
{
    int n = (int)mesh->vertices[0].feature.size();
    for (int i = 0; i < n; i++) {
        // compute max
        double max = 0.0;
        for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
            max = std::max(max, v->feature(i));
        }
        
        // normalize
        for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
            v->feature(i) /= max;
        }
    }
}

void Descriptor::compute(int n, int descriptorName)
{
    // compute descriptor
    if (descriptorName == HKS) computeHks(n);
    else if (descriptorName == WKS) computeWks(n);
    
    // normalize
    normalize();
}
