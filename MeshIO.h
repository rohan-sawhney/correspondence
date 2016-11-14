#ifndef MESH_IO_H
#define MESH_IO_H

#include <fstream>
#include "Types.h"

class MeshData;

class MeshIO {
public:
    // reads data from obj file
    static bool read(std::ifstream& in, Mesh& mesh);
    
    // reads eigenvectors and eigenvalues from file
    static void readEig(std::ifstream& in, Eigen::VectorXd& evals, Eigen::MatrixXd& evecs);
    
    // writes data in obj format
    static void write(std::ofstream& out, const Mesh& mesh);
    
    // writes eigenvectors and eigenvalues to file
    static void writeEig(std::ofstream& out, Eigen::VectorXd& evals, Eigen::MatrixXd& evecs);
    
private:
    // reserves spave for mesh vertices, uvs, normals and faces
    static void preallocateMeshElements(const MeshData& data, Mesh& mesh);
    
    // sets index for elements
    static void indexElements(Mesh& mesh);
    
    // checks if any vertex is not contained in a face
    static void checkIsolatedVertices(const Mesh& mesh);
    
    // checks if a vertex is non-manifold
    static void checkNonManifoldVertices(const Mesh& mesh);
    
    // builds the halfedge mesh
    static bool buildMesh(const MeshData& data, Mesh& mesh);
};

#endif
