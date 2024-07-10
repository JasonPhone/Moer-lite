#pragma once
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Ray/Ray.h>
// Axis-aligned Bounding box，即AABB
class AABB {
 public:
  AABB() {
    pMin = Point3f{FLT_MAX, FLT_MAX, FLT_MAX};
    pMax = Point3f{-FLT_MAX, -FLT_MAX, -FLT_MAX};
  }

  AABB(const Point3f &_pMin, const Point3f &_pMax) : pMin(_pMin), pMax(_pMax) {}

  AABB Union(const AABB &other) const;

  void Expand(const AABB &other);

  AABB Union(const Point3f &other) const;

  void Expand(const Point3f &other);

  bool Overlap(const AABB &other) const;

  bool RayIntersect(const Ray &ray, float *tMin = nullptr,
                    float *tMax = nullptr) const;

  Point3f Center() const;
  /**
   * @param p A point.
   * @return true If p is in this AABB,
   *         which means pMin <= p < pMax for each axis.
   */
  bool Contains(const Point3f &p) const;

 public:
  Point3f pMin, pMax;
};