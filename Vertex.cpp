#include "Vertex.h"
#include "HalfEdge.h"
#include "Face.h"
#include <queue>

std::vector<HalfEdge> isolated;

bool Vertex::isIsolated() const
{
    return he == isolated.begin();
}

bool Vertex::isBoundary() const
{
    HalfEdgeCIter h = he;
    do {
        if (h->onBoundary) return true;
        
        h = h->flip->next;
    } while (h != he);
    
    return false;
}

int Vertex::degree() const
{
    int k = 0;
    HalfEdgeCIter h = he;
    do {
        k++;
        
        h = h->flip->next;
    } while (h != he);
    
    return k;
}

Eigen::Vector3d Vertex::normal() const
{
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    if (isIsolated()) return normal;
    
    HalfEdgeCIter h = he;
    do {
        Eigen::Vector3d e1 = h->next->vertex->position - position;
        Eigen::Vector3d e2 = h->next->next->vertex->position - position;
        
        double d = e1.dot(e2) / sqrt(e1.squaredNorm() * e2.squaredNorm());
        if (d < -1.0) d = -1.0;
        else if (d >  1.0) d = 1.0;
        double angle = acos(d);
        
        Eigen::Vector3d n = h->face->normal();
        normal += angle * n;
        
        h = h->flip->next;
    } while (h != he);
    
    if (!normal.isZero()) normal.normalize();
    return normal;
}

double Vertex::dualArea() const
{
    double area = 0.0;
    
    HalfEdgeCIter h = he;
    do {
        area += h->face->area();
        h = h->flip->next;
        
    } while (h != he);
    
    return area / 3.0;
}

bool Vertex::isFeature(int t, int depth) const
{
    std::queue<const Vertex *> queue;
    std::unordered_map<int, bool> visited;
    
    // enqueue
    queue.push(this);
    queue.push(NULL);
    visited[index] = true;
    int levels = 0;
    
    // perform bfs
    while (!queue.empty()) {
        const Vertex *v = queue.front();
        queue.pop();
        
        if (v == NULL) {
            levels++;
            queue.push(NULL);
            if (queue.front() == NULL || levels == depth) break;
            
        } else {
            HalfEdgeCIter h = v->he;
            do {
                const Vertex *vn = &(*h->flip->vertex);
                if (!visited[vn->index]) {
                    // check if descriptor value for a particular t is greater than the neighbor's value
                    if (descriptor(t) < vn->descriptor(t)) return false;
                    
                    queue.push(vn);
                    visited[vn->index] = true;
                }
                
                h = h->flip->next;
            } while (h != v->he);
        }
    }
    
    return true;
}
