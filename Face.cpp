#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"

bool Face::isBoundary() const
{
    return he->onBoundary;
}

Eigen::Vector3d Face::normal() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    return (b-a).cross(c-a);
}

double Face::area() const
{
    if (isBoundary()) {
        return 0;
    }
    
    return 0.5 * normal().norm();
}

Eigen::Vector3d Face::centroid() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    return (a + b + c) / 3.0;
}

BoundingBox Face::boundingBox() const
{
    if (isBoundary()) {
        return BoundingBox(he->vertex->position, he->next->vertex->position);
    }
    
    BoundingBox bbox(he->vertex->position);
    bbox.expandToInclude(he->next->vertex->position);
    bbox.expandToInclude(he->next->next->vertex->position);
    
    return bbox;
}

double Face::intersect(const Eigen::Vector3d& o, const Eigen::Vector3d& d) const
{
    // Möller–Trumbore intersection algorithm
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    Eigen::Vector3d e1 = b - a;
    Eigen::Vector3d e2 = c - a;
    Eigen::Vector3d n = d.cross(e2);
    double det = e1.dot(n);
    
    // ray does not lie in the plane
    if (std::abs(det) < 1e-6) {
        return INFINITY;
    }
    
    double invDet = 1.0 / det;
    Eigen::Vector3d t = o - a;
    double u = t.dot(n) * invDet;
    
    // ray lies outside triangle
    if (u < 0.0 || u > 1.0) {
        return INFINITY;
    }
    
    Eigen::Vector3d q = t.cross(e1);
    double v = d.dot(q) * invDet;
    
    // ray lies outside the triangle
    if (v < 0.0 || v + u > 1.0) {
        return INFINITY;
    }
    
    // check for intersection
    double s = e2.dot(q) * invDet;
    if (s > 0) {
        return s;
    }
    
    // no hit
    return INFINITY;
}
