#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Task 1.
    int n = points.size();

    if(n == 1){
      return points;

    }else{
      std::vector<Vector2D> next(n-1);

      for(int i = 0; i < n-1; ++i){
        next[i] = (1-t) * points[i] + t * points[i+1];
      }

      return next;
    }

  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Task 2.

    int n = points.size();

    if(n == 1){
      return points;

    }else{
      std::vector<Vector3D> next(n-1);

      for(int i = 0; i < n-1; ++i){
        next[i] = (1-t) * points[i] + t * points[i+1];
      }

      return next;
    }

  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Task 2.
    int n = points.size();
    std::vector<Vector3D> step_results = points;

    while(n > 1){
      step_results = evaluateStep(step_results, t);
      n--;
    }

    return step_results[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Task 2.
    int n = controlPoints.size();
    std::vector<Vector3D> points(n);

    for(int i = 0; i < n; ++i){
      points[i] = evaluate1D(controlPoints[i], u);
    }

    return evaluate1D(points, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Task 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    HalfedgeCIter h = halfedge();
    VertexCIter v = h->vertex();
    Vector3D center = v->position;
    std::vector<Vector3D> neighbors;

    do {
        HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
        VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                                          // h->vertex() is v, whereas h_twin->vertex()
                                          // is the neighboring vertex
        neighbors.push_back(v->position);
        h = h_twin->next();               // move to the next outgoing half-edge of the vertex
    } while(h != halfedge());

    int n = neighbors.size();
    std::vector<double> areas(n);
    double area_sum = 0;
    std::vector<Vector3D> normals(n);
    for(int i = 0; i < n; ++i){
      Vector3D p0 = neighbors[i];
      Vector3D p1 = neighbors[(i+1)%n];

      double edge0 = (p0 - p1).norm();
      double edge1 = (p1 - center).norm();
      double edge2 = (center - p0).norm();

      double p = (edge0 + edge1 + edge2) / 2;
      double area = sqrt(p * (p-edge0) * (p-edge1) * (p-edge2));

      areas[i] = area;
      area_sum += area;

      Vector3D v0 = p0 - center;
      Vector3D v1 = p1 - center;
      Vector3D normal = cross(v1, v0);
      normal.normalize();

      normals[i] = normal;
    }

    Vector3D N;
    for(int i = 0; i < n; ++i){
      N += areas[i] / area_sum * normals[i];
    }

    return N;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Task 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    if(!e0->isBoundary()){

      // COLLECT ELEMENTS
      
      // half edges
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      // vertices
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();

      // edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      // faces
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();


      // REASSIGN ELEMENTS

      // half edges
      // inside
      h0->setNeighbors(h1, h3, v3, e0, f0);
      h1->setNeighbors(h2, h7, v2, e2, f0);
      h2->setNeighbors(h0, h8, v0, e3, f0);
      h3->setNeighbors(h4, h0, v2, e0, f1);
      h4->setNeighbors(h5, h9, v3, e4, f1);
      h5->setNeighbors(h3, h6, v1, e1, f1);
      // outside
      // h6->setNeighbors(h6->next(), h5, v2, e1, h6->face());
      // h7->setNeighbors(h7->next(), h1, v0, e2, h7->face());
      // h8->setNeighbors(h8->next(), h2, v3, e3, h8->face());
      // h9->setNeighbors(h9->next(), h4, v1, e4, h9->face());
      h6->twin() = h5;
      h7->twin() = h1;
      h8->twin() = h2;
      h9->twin() = h4;

      // vertices
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h3;
      v3->halfedge() = h0;

      // edges
      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h1;
      e3->halfedge() = h2;
      e4->halfedge() = h4;

      // faces
      f0->halfedge() = h0;
      f1->halfedge() = h3;
    }
    
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Task 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    VertexIter v4 = newVertex(); // the output vertex iter

    if(!e0->isBoundary()){ // non boundary case

      // COLLECT ELEMENTS
      
      // half edges
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      // vertices
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();

      // edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();

      // faces
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();


      // ALLOCATE NEW ELEMENTS

      // half edges
      HalfedgeIter h10 = newHalfedge();
      HalfedgeIter h11 = newHalfedge();
      HalfedgeIter h12 = newHalfedge();
      HalfedgeIter h13 = newHalfedge();
      HalfedgeIter h14 = newHalfedge();
      HalfedgeIter h15 = newHalfedge();

      // vertices
      // VertexIter v4 = newVertex();
      v4->position = (v0->position + v1->position) / 2;

      // edges
      EdgeIter e5 = newEdge();
      EdgeIter e6 = newEdge();
      EdgeIter e7 = newEdge();

      // faces
      FaceIter f2 = newFace();
      FaceIter f3 = newFace();


      // REASSIGN ELEMENTS

      // half edges
      // inside
      h0->setNeighbors(h1, h14, v0, e0, f0);
      h1->setNeighbors(h2, h3, v4, e5, f0);
      h2->setNeighbors(h0, h7, v2, e2, f0);

      h3->setNeighbors(h4, h1, v2, e5, f1);
      h4->setNeighbors(h5, h10, v4, e6, f1);
      h5->setNeighbors(h3, h6, v1, e1, f1);

      h10->setNeighbors(h11, h4, v1, e6, f2);
      h11->setNeighbors(h12, h13, v4, e7, f2);
      h12->setNeighbors(h10, h9, v3, e4, f2);

      h13->setNeighbors(h14, h11, v3, e7, f3);
      h14->setNeighbors(h15, h0, v4, e0, f3);
      h15->setNeighbors(h13, h8, v0, e3, f3);

      // outside
      h6->twin() = h5;
      h7->twin() = h2;
      h8->twin() = h15;
      h9->twin() = h12;

      // vertices
      v0->halfedge() = h0;
      v1->halfedge() = h10;
      v2->halfedge() = h3;
      v3->halfedge() = h13;
      v4->halfedge() = h1;

      // edges
      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h2;
      e3->halfedge() = h15;
      e4->halfedge() = h12;
      e5->halfedge() = h1;
      e6->halfedge() = h4;
      e7->halfedge() = h11;

      // faces
      f0->halfedge() = h0;
      f1->halfedge() = h3;
      f2->halfedge() = h10;
      f3->halfedge() = h13;
    
    }else{ // boundary case

      // COLLECT ELEMENTS

      // half edges
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();

      // vertices
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h1->vertex();
      VertexIter v2 = h2->vertex();

      // edges
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();

      // faces
      FaceIter f0 = h0->face();


      // ALLOCATE NEW ELEMENTS

      // half edges
      HalfedgeIter h6 = newHalfedge();
      HalfedgeIter h7 = newHalfedge();
      HalfedgeIter h8 = newHalfedge();
      HalfedgeIter h9 = newHalfedge();

      // vertices
      // VertexIter v4 = newVertex();
      v4->position = (v0->position + v1->position) / 2;

      // edges
      EdgeIter e3 = newEdge();
      EdgeIter e4 = newEdge();

      // faces
      FaceIter f1 = newFace();


      // REASSIGN ELEMENTS

      // half edges
      // inside
      h0->setNeighbors(h1, h6, v0, e0, f0);
      h1->setNeighbors(h2, h7, v4, e4, f0);
      h2->setNeighbors(h0, h4, v2, e2, f0);

      h7->setNeighbors(h8, h1, v2, e4, f1);
      h8->setNeighbors(h9, h3, v4, e3, f1);
      h9->setNeighbors(h7, h5, v1, e1, f1);

      // outside
      h4->twin() = h2;
      h5->twin() = h9;

      h6->setNeighbors(h3->next(), h0, v4, e0, h3->face());
      h3->setNeighbors(h6, h8, v1, e3, h3->face());
      
      // vertices
      v0->halfedge() = h0;
      v1->halfedge() = h9;
      v2->halfedge() = h2;
      v4->halfedge() = h1;

      // edges
      e0->halfedge() = h0;
      e1->halfedge() = h9;
      e2->halfedge() = h2;
      e3->halfedge() = h8;
      e4->halfedge() = h1;

      // faces
      f0->halfedge() = h0;
      f1->halfedge() = h7;
      
    }

    return v4;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Task 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    for( VertexIter v1 = mesh.verticesBegin(); v1 != mesh.verticesEnd(); v1++ ){
      HalfedgeCIter h = v1->halfedge();
      Vector3D neighbor_pos_sum(0, 0, 0);
      int n = 0;
      do {
        HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
        VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
                                          // h->vertex() is v, whereas h_twin->vertex()
                                          // is the neighboring vertex
        neighbor_pos_sum += v->position;
        n++;
        h = h_twin->next();               // move to the next outgoing half-edge of the vertex
      } while(h != v1->halfedge());
      double u = ((n == 3) ? 0.1875 : 0.375/n);
      v1->newPosition = (1-n*u) * v1->position + u * neighbor_pos_sum;
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    for( EdgeIter e2 = mesh.edgesBegin(); e2 != mesh.edgesEnd(); e2++ ){
      if(!e2->isBoundary()){
        HalfedgeCIter h = e2->halfedge();
        VertexCIter A = h->vertex();
        VertexCIter B = h->twin()->vertex();
        VertexCIter C = h->next()->next()->vertex();
        VertexCIter D = h->twin()->next()->next()->vertex();
        e2->newPosition = (A->position + B->position) * 3/8 + (C->position + D->position) / 8;
      }else{
        HalfedgeCIter h = e2->halfedge();
        if(h->face()->isBoundary()){
          h = h->twin();
        }
        VertexCIter A = h->vertex();
        VertexCIter B = h->next()->vertex();
        e2->newPosition = (A->position + B->position) / 2;
      }
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    for( EdgeIter e3 = mesh.edgesBegin(); e3 != mesh.edgesEnd(); e3++ ){
      if(!e3->isNew){
        HalfedgeCIter hc = e3->halfedge();
        if(!hc->vertex()->isNew && !hc->twin()->vertex()->isNew){
          VertexIter v = mesh.splitEdge(e3);
          v->isNew = true;
          v->newPosition = e3->newPosition;
          // this part depends on the design of the edge split function !!!
          HalfedgeIter h = v->halfedge();
          if(!e3->isBoundary()){
            h->edge()->isNew = true;
            h->twin()->next()->twin()->next()->edge()->isNew = true;
          }else{
            h->edge()->isNew = true;
          }
        }
      }
    }
    
    // 4. Flip any new edge that connects an old and new vertex.
    for( EdgeIter e4 = mesh.edgesBegin(); e4 != mesh.edgesEnd(); e4++ ){
      if(e4->isNew){
        HalfedgeCIter h = e4->halfedge();
        VertexCIter v0 = h->vertex();
        VertexCIter v1 = h->twin()->vertex();
        if(v0->isNew ^ v1->isNew == 1) mesh.flipEdge(e4);
      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for( VertexIter v5 = mesh.verticesBegin(); v5 != mesh.verticesEnd(); v5++ ){
      v5->position = v5->newPosition;
      v5->isNew = false;
    }
    // reset the isNew flags
    for( EdgeIter e5 = mesh.edgesBegin(); e5 != mesh.edgesEnd(); e5++ ){
      e5->isNew = false;
    }
  }
}
