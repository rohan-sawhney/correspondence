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

void Descriptor::computeEig(const Eigen::SparseMatrix<double>& W,
                            const Eigen::SparseMatrix<double>& A)
{
    Spectra::SparseSymMatProd<double> opW(W);
    Spectra::SparseCholesky<double> opA(A);
    
    Spectra::SymGEigsSolver<double,
                            Spectra::SMALLEST_MAGN,
                            Spectra::SparseSymMatProd<double>,
                            Spectra::SparseCholesky<double>,
                            Spectra::GEIGS_CHOLESKY> geigs(&opW, &opA, K, 2*K);
    
    geigs.init();
    geigs.compute();
    
    if (geigs.info() == Spectra::SUCCESSFUL) {
        evals = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    
    } else {
        std::cout << "Eigen computation failed" << std::endl;
    }
    
    Eigen::MatrixXd err = W*evecs - A*evecs*evals.asDiagonal();
    std::cout << "||Lx - λAx||_inf = " << err.array().abs().maxCoeff() << std::endl;
}

void buildAdjacency(Mesh *mesh, Eigen::SparseMatrix<double>& W)
{
    std::vector<Eigen::Triplet<double>> WTriplets;
    
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        do {
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan());
            sumCoefficients += coefficient;
            
            WTriplets.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -coefficient));
            
            he = he->flip->next;
        } while (he != v->he);
        
        WTriplets.push_back(Eigen::Triplet<double>(v->index, v->index, sumCoefficients));
    }
    
    W.setFromTriplets(WTriplets.begin(), WTriplets.end());
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
        
        // build adjacency operator
        Eigen::SparseMatrix<double> W(v, v);
        buildAdjacency(mesh, W);
        
        // build area matrix
        Eigen::SparseMatrix<double> A(v, v);
        buildAreaMatrix(mesh, A);
        
        // compute eigenvectors and eigenvalues
        computeEig(W, A);
        
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

void extrapolateEvals(double& xhat, double& yhat, double& m, const Eigen::VectorXd& evals)
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

double computeMaxEigenvalue(const Eigen::SparseMatrix<double>& L)
{
    Spectra::SparseSymMatProd<double> opL(L);
    
    Spectra::SymEigsSolver<double,
    Spectra::LARGEST_MAGN,
    Spectra::SparseSymMatProd<double>> eigs(&opL, 1, 3);
    
    eigs.init();
    eigs.compute();
    
    if (eigs.info() == Spectra::SUCCESSFUL) return eigs.eigenvalues()[0];
    else std::cout << "Eigen computation failed" << std::endl;
    
    return 0.0;
}

void computeLaplacians(std::vector<Eigen::SparseMatrix<double>>& Ls,
                       std::vector<Eigen::SparseMatrix<double>>& As,
                       std::vector<Eigen::SparseMatrix<double>>& Ainvs,
                       std::vector<double>& maxEvals,
                       std::vector<double>& conditionNumbersA,
                       const MultiresMesh& mrm)
{
    int lods = mrm.numLods();
    
    for (int l = 0; l < lods; l++) {
        Mesh *lod = mrm.lod(l);
        int v = (int)lod->vertices.size();
        double mina = INFINITY, maxa = -INFINITY;
        
        // resize
        Ls[l].resize(v, v); As[l].resize(v, v); Ainvs[l].resize(v, v);
        
        // compute laplacian and area matrices
        buildAdjacency(lod, Ls[l]);
        buildAreaMatrix(lod, As[l]);
        for (VertexIter v = lod->vertices.begin(); v != lod->vertices.end(); v++) {
            double a = As[l].coeffRef(v->index, v->index);
            Ainvs[l].insert(v->index, v->index) = 1.0/a;
            maxa = std::max(maxa, a);
            mina = std::min(mina, a);
        }
        Ainvs[l].makeCompressed();
        Ls[l] = Ainvs[l]*Ls[l];
        
        // compute max evals and condition numbers
        maxEvals[l] = computeMaxEigenvalue(Ls[l]);
        conditionNumbersA[l] = sqrt(maxa / mina);
    }
}

double computeBinomialElementNorm(const double lambda, const double a, const double t, const int m)
{
    double binomial = 0.0;
    double f = 1.0;
    for (int k = 0; k <= m-1; k++) {
        binomial *= (lambda - k);
        f *= (k+1);
    }
    binomial /= f;
    
    return fabs(binomial*(exp(-t) - 1))*a;
}

void computeBinomialRepresentation(Eigen::SparseMatrix<double>& Kt,
                                   const Eigen::SparseMatrix<double>& L,
                                   const double t, const int iters)
{
    int k = 0;
    double mf = 1;
    Kt.setZero();
    Eigen::SparseMatrix<double> Id(L.cols(), L.cols()); Id.setIdentity();
    Eigen::SparseMatrix<double> Q = Id;
    for (int m = 0; m <= iters; m++) {
        if (k <= m-1) {
            Q = Q*(L - k*Id);
            k++;
        }
        
        if (m > 0) mf *= m;
        Q /= mf;
        Kt += Q*pow(exp(-t) - 1, m);
    }
}

void sparsify(Eigen::SparseMatrix<double>& Kt, double eps)
{
    for (int i = 0; i < Kt.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Kt, i); it; ++it) {
            if (it.valueRef() < eps) it.valueRef() = 0.0;
        }
    }
    Kt.prune(0.0);
}

void Descriptor::computeFastHks()
{
    const int d = 2;
    const int C = 1000;
    const double c = 0.2;
    const double reps = 1e-4;
    const double seps = 1e-6;
    const int iters = 15;
    const std::vector<int> ts = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    
    // extrapolate eigenvalues
    double xhat, yhat, m;
    extrapolateEvals(xhat, yhat, m, evals);
    
    // 1. build mutliresolution structure
    MultiresMesh mrm(mesh, d, C);
    mrm.build();
    
    // build laplacian and area matrices and compute max evals for laplacians and condition numbers for A
    int lods = mrm.numLods();
    std::vector<Eigen::SparseMatrix<double>> Ls(lods), As(lods), Ainvs(lods);
    std::vector<double> maxEvals(lods), conditionNumbersA(lods);
    computeLaplacians(Ls, As, Ainvs, maxEvals, conditionNumbersA, mrm);
    
    for (int i = 0; i < N; i++) {
        // 2. choose coarsest resolution level l s.t. cn > r(t)
        int r = 0, l = lods-1;
        while (exp(-(yhat - m*(r - xhat))*ts[i]) > reps) r++;
        while (l > 0 && c*mrm.lod(l)->vertices.size() < r) l--;
        
        // 3. compute sparse heat kernel on l
        int v = (int)mrm.lod(l)->vertices.size();
        Eigen::SparseMatrix<double> Kt(v, v);
        double t = 0.0, t1 = 0.05; // NOTE: I'm not sure if t1 is being set correctly 
        int s = 0;
        
        // compute Kt for small t
        while (computeBinomialElementNorm(maxEvals[l], conditionNumbersA[l], t1, iters) > seps) t1 += 0.05;
        while ((t = ts[i]/pow(2, s)) > t1) s++;
        computeBinomialRepresentation(Kt, Ls[l], t, iters);

        // compute Kt for ts[i]
        for (int j = 0; j < s; j++) {
            sparsify(Kt, seps);
            Eigen::SparseMatrix<double> KtA = Kt*As[l];
            Kt = KtA*KtA*Ainvs[l];
            t = 2*t;
        }
        
        // 4. project sparse heat kernel on finest resolution level
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
