#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bbox;
  Vector3D center_sum; // to calculate the average
  size_t size = 0; // number of Primitives in this node

  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
    center_sum += bb.centroid();
    size++;
  }

  BVHNode *node = new BVHNode(bbox);
  node->start = start;
  node->end = end;

  if(size > max_leaf_size){
    // find the longest axis to split
    int axis = (bbox.extent.x > bbox.extent.y) ? 0 : 1;
    double extent = (axis == 0) ? bbox.extent.x : bbox.extent.y;
    axis = (extent > bbox.extent.z) ? axis : 2;

    std::vector<Primitive *>* left = new std::vector<Primitive *>;
    std::vector<Primitive *>* right = new std::vector<Primitive *>;

    // split into left & right children
    for (auto p = start; p != end; p++) {
      Vector3D c = (*p)->get_bbox().centroid();
      if(c[axis] < center_sum[axis]/size && left->size() < size * 0.7){ // less than average && not half size yet
        left->push_back(*p);
      }else if(c[axis] > center_sum[axis]/size && right->size() < size * 0.7){
        right->push_back(*p);
      }else{
        if(left->size() < right->size()){
          left->push_back(*p);
        }else{
          right->push_back(*p);
        }
      }
    }

    node->l = construct_bvh(left->begin(), left->end(), max_leaf_size);
    node->r = construct_bvh(right->begin(), right->end(), max_leaf_size);
  }

  return node;

}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  double t0, t1;
  total_isects++;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (t0 > ray.max_t || t1 < ray.min_t) {
    return false;
  }

  if(node->isLeaf()){
    for (auto p : primitives) {
      total_isects++;
      if (p->has_intersection(ray))
        return true;
    }
    return false;
  }

  return has_intersection(ray, node->l) || has_intersection(ray, node->r);

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.

  double t0, t1;
  total_isects++;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (t0 > ray.max_t || t1 < ray.min_t) {
    return false;
  }
    
  if(node->isLeaf()){
    bool hit = false;
    for (auto p : primitives) {
      total_isects++;
      hit = p->intersect(ray, i) || hit;
    }
    return hit;
  }

  return intersect(ray, i, node->l) || intersect(ray, i, node->r);

}

} // namespace SceneObjects
} // namespace CGL
