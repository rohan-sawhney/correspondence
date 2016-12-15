#include "Mesh.h"
#include "MeshIO.h"

Mesh::Mesh()
{
    
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        name = fileName;
        normalize();
    }
    
    in.close();
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    out.close();
    return false;
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin and determine radius
    double rMax = 0;
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}


double Mesh::computeGaussCurvature(Eigen::VectorXd& K)
{
    double maxGauss = -INFINITY;
    for (size_t i = 0; i < vertices.size(); i++) {
        K(i) = vertices[i].angleDefect() / vertices[i].dualArea();
        
        if (maxGauss < fabs(K(i))) maxGauss = fabs(K(i));
    }
    
    return maxGauss;
}

void Mesh::buildLaplacian(Eigen::SparseMatrix<double>& L) const
{
    std::vector<Eigen::Triplet<double>> LTriplet;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double dualArea = v->dualArea();
        double sumCoefficients = 0.0;
        do {
            // (cotA + cotB) / 2A
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan()) / dualArea;
            sumCoefficients += coefficient;
            
            LTriplet.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -coefficient));
            
            he = he->flip->next;
        } while (he != v->he);

        LTriplet.push_back(Eigen::Triplet<double>(v->index, v->index, sumCoefficients));
    }
    
    L.setFromTriplets(LTriplet.begin(), LTriplet.end());
}

void Mesh::buildSimpleAverager(Eigen::SparseMatrix<double>& L) const
{
    std::vector<Eigen::Triplet<double>> LTriplet;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double degree = v->degree();
        do {
            LTriplet.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, 1.0/degree));
            he = he->flip->next;
        } while (he != v->he);
    }
    
    L.setFromTriplets(LTriplet.begin(), LTriplet.end());
}


double Mesh::computeMeanCurvature(Eigen::VectorXd& H)
{
    Eigen::SparseMatrix<double> L((int)vertices.size(), (int)vertices.size());
    buildLaplacian(L);
    
    Eigen::MatrixXd x;
    x.resize((int)vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++) {
        x.row(i) = vertices[i].position;
    }
    x = L * x;

    // set absolute mean curvature
    double maxMean = -INFINITY;
    for (size_t i = 0; i < vertices.size(); i++) {
        H(i) = 0.5 * x.row(i).norm();
        
        if (maxMean < H(i)) maxMean = H(i);
    }
    
    return maxMean;
}

void Mesh::computeCurvatures()
{
    int v = (int)vertices.size();
    
    Eigen::VectorXd K(v);
    double maxGauss = computeGaussCurvature(K);
    
    Eigen::VectorXd H(v);
    double maxMean = computeMeanCurvature(H);
    
    // compute principal curvatures and normalize gauss, mean curvature
    for (int i = 0; i < v; i++) {
        double dis = sqrt(H(i)*H(i) - K(i));
        
        vertices[i].k1 = H(i) + dis;
        vertices[i].k2 = H(i) - dis;
        
        vertices[i].gaussCurvature = K(i) / maxGauss;
        vertices[i].meanCurvature = H(i) / maxMean;
    }
}

