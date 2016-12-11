#ifndef MESH_H
#define MESH_H

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"

class Mesh {
public:
    // default constructor
    Mesh();
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
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
};

#endif
