#ifndef MESH_H
#define MESH_H

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <Eigen/SparseCore>

class Mesh {
public:
    // default constructor
    Mesh();
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
   
    // computes principal, gaussian and mean curvatures
    void computeCurvatures();
 
    void buildSimpleAverager(Eigen::SparseMatrix<double>& L) const;

    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;
    std::string name;

private:
    // center mesh about origin and rescale to unit radius
    void normalize();

    // computes gaussian curvature per vertex
    double computeGaussCurvature(Eigen::VectorXd& K);
    
    // builds Laplace Beltrami operator
    void buildLaplacian(Eigen::SparseMatrix<double>& L) const;
    
    // computes mean curvature per vertex
    double computeMeanCurvature(Eigen::VectorXd& H);
};

#endif
