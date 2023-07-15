#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.
  
  // Moller Trumbore Algr.
  Vector3D S, S1, S2, E1, E2;
  E1 = p2 - p1; 
  E2 = p3 - p1; 
  S = r.o - p1; 
  S1 = cross(r.d, E2);
  S2 = cross(S, E1);
  double t = dot(S2, E2) / dot(S1, E1);
  double b2 = dot(S1, S) / dot(S1, E1);
  double b3 = dot(S2, r.d) / dot(S1, E1);

  if(t > r.min_t && t < r.max_t && b2 > 0 && b3 > 0 && (1-b2-b3 > 0)) return true;

  return false;

}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

  Vector3D S, S1, S2, E1, E2;
  E1 = p2 - p1; 
  E2 = p3 - p1; 
  S = r.o - p1; 
  S1 = cross(r.d, E2);
  S2 = cross(S, E1);
  double t = dot(S2, E2) / dot(S1, E1);
  double b2 = dot(S1, S) / dot(S1, E1);
  double b3 = dot(S2, r.d) / dot(S1, E1);
  double b1 = 1 - b2 - b3;

  if(t > r.min_t && t < r.max_t && b1 > 0 && b2 > 0 && b3 > 0){
    r.max_t = t;
    isect->t = t;
    isect->n = b1 * n1 + b2 * n2 + b3 * n3;
    isect->primitive = this;
    isect->bsdf = this->get_bsdf();
    return true;
  }

  return false;
}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
