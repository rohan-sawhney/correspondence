#ifndef MULTI_RES_MESH_H
#define MULTI_RES_MESH_H

#include "Types.h"
#define _USE_MATH_DEFINES
#include <openmesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#undef _USE_MATH_DEFINES

struct Traits : public OpenMesh::DefaultTraits {
    // double valued point
    typedef OpenMesh::Vec3d Point;
    
    // not storing previous halfedge
    HalfedgeAttributes(OpenMesh::Attributes::None);
};

typedef OpenMesh::TriMesh_ArrayKernelT<Traits> OMesh;

class MultiresMesh {
public:
    // construction
    MultiresMesh(Mesh *mesh0, double d0, double C0);
    
    // destructor
    ~MultiresMesh();
    
    // build
    void build();
    
    // returns lod count
    int numLods();
    
    // returns l'th lod
    Mesh* lod(int l);
    
private:
    // loads mesh into oMesh
    void buildOMesh(Mesh *mesh);
    
    // builds mesh from decimated oMesh
    void buildMesh(Mesh *mesh);
    
    // decimates
    Mesh* decimate(int v);
    
    // projects mesh1 onto mesh2
    void project(Mesh *mesh1, Mesh *mesh2);
    
    // Member variables
    std::vector<Mesh *> lods;
    OMesh oMesh;
    double d;
    double C;
};

#endif
