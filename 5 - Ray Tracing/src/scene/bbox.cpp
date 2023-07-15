#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

void swap(double& d1, double& d2) {
	double temp = d1;
	d1 = d2;
	d2 = temp;
}

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double min_xt = (min.x - r.o.x) / r.d.x;
  double max_xt = (max.x - r.o.x) / r.d.x;
  if (r.d.x < 0) swap(min_xt, max_xt);

  double min_yt = (min.y - r.o.y) / r.d.y;
  double max_yt = (max.y - r.o.y) / r.d.y;
  if (r.d.y < 0) swap(min_yt, max_yt);

  double min_zt = (min.z - r.o.z) / r.d.z;
  double max_zt = (max.z - r.o.z) / r.d.z;
  if (r.d.z < 0) swap(min_zt, max_zt);

  double min_t = std::max(std::max(min_xt, min_yt), min_zt);
  double max_t = std::min(std::min(max_xt, max_yt), max_zt);
  
  if (min_t < max_t) {
    t0 = min_t;
    t1 = max_t;
    return true;
  }
  return false;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
