#include "Descriptor.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "MultiresMesh.h"
#include <spectra/include/MatOp/SparseSymMatProd.h>
#include <spectra/include/MatOp/SparseCholesky.h>
#include <spectra/include/SymGEigsSolver.h>
#define N 10

Descriptor::Descriptor(Mesh *mesh0):
mesh(mesh0)
{
    
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
        
        WTriplets.push_back(Eigen::Triplet<double>(v->index, v->index,
                                                   sumCoefficients + 1e-8));
    }
    
    W.setFromTriplets(WTriplets.begin(), WTriplets.end());
}

void buildAreaMatrix(Mesh *mesh, Eigen::SparseMatrix<double>& A, const double scale)
{
    std::vector<Eigen::Triplet<double>> ATriplets;
    
    double sum = 0.0;
    for (VertexCIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        double area = v->dualArea();
        ATriplets.push_back(Eigen::Triplet<double>(v->index, v->index, area));
        sum += area;
    }
    
    A.setFromTriplets(ATriplets.begin(), ATriplets.end());
    A *= scale/sum;
}

void Descriptor::computeEig(const int K, const bool loadEig)
{
    // read eigvalues and eigenvectors from file
    std::string filename = mesh->name;
    filename.replace(filename.find_last_of(".")+1, 3, "eig");
    std::ifstream in(filename);
    
    if (loadEig && in.is_open()) {
        MeshIO::readEig(in, evals, evecs);
        in.close();
        
    } else {
        int v = (int)mesh->vertices.size();
        
        // build adjacency operator
        Eigen::SparseMatrix<double> W(v, v);
        buildAdjacency(mesh, W);
        
        // build area matrix
        Eigen::SparseMatrix<double> A(v, v);
        buildAreaMatrix(mesh, A, mesh->vertices.size());
        
        // compute eigenvectors and eigenvalues
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
        
        // write eigvalues and eigenvectors to file
        std::ofstream out(filename);
        if (out.is_open()) {
            MeshIO::writeEig(out, evals, evecs);
            out.close();
        }
    }
}

void Descriptor::computeHks()
{
    const int K = (int)evals.size();
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
    const int K = (int)evals.size();
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

void computeLaplacians(std::vector<Eigen::SparseMatrix<double>>& Ls,
                       std::vector<Eigen::SparseMatrix<double>>& As,
                       const MultiresMesh& mrm)
{
    double scale = mrm.lod(0)->vertices.size();
    for (int l = 0; l < mrm.numLods(); l++) {
        Mesh *lod = mrm.lod(l);
        int v = (int)lod->vertices.size();
        
        // resize
        Ls[l].resize(v, v); As[l].resize(v, v);
        
        // compute laplacian and area matrices
        buildAdjacency(lod, Ls[l]);
        buildAreaMatrix(lod, As[l], scale);
        Ls[l] = As[l].cwiseInverse()*Ls[l];
    }
}

void computeBinomialEntries(std::vector<Eigen::SparseMatrix<double>>& binomialSeries,
                            const Eigen::SparseMatrix<double>& L)
{
    int k = 0;
    Eigen::SparseMatrix<double> Id(L.cols(), L.cols()); Id.setIdentity();
    Eigen::SparseMatrix<double> Q = Id;
    for (int m = 0; m < (int)binomialSeries.size(); m++) {
        if (k <= m-1) {
            Q = Q*(L - k*Id);
            k++;
        }
        
        if (m > 0) Q /= m;
        binomialSeries[m] = Q;
    }
}

void computeExponentialRepresentation(Eigen::SparseMatrix<double>& Kt, const double t,
                                      const std::vector<Eigen::SparseMatrix<double>>& binomialSeries)
{
    Kt.setZero();
    for (int m = 0; m < (int)binomialSeries.size(); m++) {
        Kt += binomialSeries[m]*pow(exp(-t) - 1, m);
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
    const int bN = 15;
    const std::vector<int> ts = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    
    // initialize descriptors
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->descriptor = Eigen::VectorXd::Zero(N);
    }
    
    // extrapolate eigenvalues
    double xhat, yhat, m;
    extrapolateEvals(xhat, yhat, m, evals);
    
    // 1. build mutliresolution structure
    MultiresMesh mrm(mesh, d, C);
    mrm.build();
    
    // build laplacian and area matrices
    int lods = mrm.numLods();
    std::vector<Eigen::SparseMatrix<double>> Ls(lods), As(lods);
    computeLaplacians(Ls, As, mrm);
    
    for (int i = 0; i < N; i++) {
        // 2. choose coarsest resolution level l s.t. cn > r(t)
        int r = 0, l = lods-1;
        while (exp(-(yhat - m*(r - xhat))*ts[i]) > reps) r++;
        while (l > 0 && c*mrm.lod(l)->vertices.size() < r) l--;
        std::cout << "l: " << l << " r: " << r << std::endl;
        
        // 3. compute sparse heat kernel on l
        int v = (int)mrm.lod(l)->vertices.size();
        Eigen::SparseMatrix<double> Kt(v, v);
        std::vector<Eigen::SparseMatrix<double>> binomialSeries(bN+1);
        double t1 = 0.001, t = 0.0;
        int s = 0;
        
        // compute Kt for small t
        computeBinomialEntries(binomialSeries, Ls[l]);
        while ((binomialSeries[bN]*pow(exp(-t1) - 1, bN)).norm() < seps) t1 += 0.001;
        while ((t = ts[i]/pow(2, s)) > t1) s++;
        computeExponentialRepresentation(Kt, t, binomialSeries);
        std::cout << "t1: " << t1 << " t: " << t << " s: " << s << std::endl;
        
        // compute Kt for ts[i]
        for (int j = 0; j < s; j++) {
            sparsify(Kt, seps);
            Kt = Kt*Kt;
            std::cout << "Kt 1: " << Kt.toDense().lpNorm<1>() << std::endl;
            std::cout << "Kt 2: " << Kt.norm() << std::endl;
        }
        
        // 4. project sparse heat kernel on finest resolution level
        Eigen::SparseMatrix<double> P = mrm.prolongationMatrix(l);
        std::cout << "P 1: " << P.toDense().lpNorm<1>() << std::endl;
        std::cout << "P 2: " << P.norm() << std::endl;
        Kt = P*Kt*P.transpose();
        std::cout << "pKt 1: " << Kt.toDense().lpNorm<1>() << std::endl;
        std::cout << "pKt 2: " << Kt.norm() << std::endl;
        
        // set descriptor value for current time step
        for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
            v->descriptor(i) = Kt.coeffRef(v->index, v->index);
        }
    }
}

void Descriptor::computeWks()
{
    const int K = (int)evals.size();
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

void Descriptor::computeCurve() {

	std::cout << "Computing curvature descriptor..." << std::endl;
   
	std::vector<int> smoothLevels = {0, 1, 5, 20, 50};
	unsigned int dSize = 2*smoothLevels.size(); // mean,gaussian at each level
	 

	// === Preparation
	
	// Compute the actual curvature values
	mesh->computeCurvatures();

	// Allocate space for the descriptors
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
        v->descriptor = Eigen::VectorXd::Zero(dSize);
	}


	// Push the curvatures to vectors
	unsigned int nVerts = mesh->vertices.size();
    Eigen::VectorXd currGauss(nVerts);
    Eigen::VectorXd currMean(nVerts);
	size_t i = 0;
    for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
		currGauss(i) = v->gaussCurvature;
		currMean(i) = v->meanCurvature;
		i++;
	}

	// Compute a shifted-laplacian matrix for taking averages
    Eigen::SparseMatrix<double> avgM(nVerts, nVerts);
    mesh->buildSimpleAverager(avgM);

	// === Smoothing and saving

	// For each of the smoothing levels, smooth an appropriate number of times, then save the descriptor
	int smoothingStepsCompleted = 0;
	int iLevel = 0;
	for(int smoothLevel : smoothLevels) {

		// Smooth as needed
		while(smoothingStepsCompleted < smoothLevel) {

			currGauss = avgM * currGauss;
			currMean = avgM * currMean;

			smoothingStepsCompleted++;
		}


		// Save
		size_t iVert = 0;
		for (VertexIter v = mesh->vertices.begin(); v != mesh->vertices.end(); v++) {
			v->descriptor(iLevel + 0) = currGauss(iVert);
			v->descriptor(iLevel + 1) = currMean(iVert);
			iVert++;
		}
		
		iLevel += 2;
	}

	std::cout << "... Done computing curvature descriptor." << std::endl;
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

void Descriptor::compute(int descriptor, bool loadEig, std::string outFilename)
{
    // compute descriptor
    std::string descriptorName;
    switch (descriptor) {
        case HKS:
            computeEig(300, loadEig);
            computeHks();
            descriptorName = "hks";
            break;
        case FAST_HKS:
            computeEig(30, loadEig);
            computeFastHks();
            descriptorName = "fks";
            break;
        case WKS:
            computeEig(300, loadEig);
            computeWks();
            descriptorName = "wks";
            break;
        case CURVE:
            computeCurve();
            descriptorName = "curve";
            break;
    }
    
    // normalize
    normalize();
    
    // write to file
    //std::string filename = mesh->name;
    //filename.replace(filename.find_last_of(".")+1, 3, descriptorName);
    std::ofstream out(outFilename);
    
    if (out.is_open()) {
        MeshIO::writeDescriptor(out, *mesh);
        out.close();
    } else {
		std::cout << "Not writing descriptor, no valid path specified" << std::endl;
	}
}
