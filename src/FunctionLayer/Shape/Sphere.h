#pragma once
#include "CoreLayer/Math/Constant.h"
#include "CoreLayer/Math/Geometry.h"
#include "FastMath/FastMath.h"
#include "Shape.h"
#include <math.h>
// TODO 当前只有Transform中的translate对sphere生效
class Sphere : public Shape {
public:
  Sphere() = delete;

  Sphere(const Json &json);

  virtual bool rayIntersectShape(Ray &ray, int *primID, float *u,
                                 float *v) const override;

  virtual void fillIntersection(float distance, int primID, float u, float v,
                                Intersection *intersection) const override;
  virtual void uniformSampleOnSurface(Vector2f sample,
                                      Intersection *intersection,
                                      float *pdf) const override {
    // TODO finish this
    /*
      c = 1 / 4pi
      pdf(t, p) = c sint
      pdf(t) = c sint * 2pi = 0.5 sint
      pdf(p) = [0, pi] c sint dt = -c [0, pi] dcost = 2c = 1 / 2pi

      cdf(t) = 0.5 [0, x] sint dt = -0.5 [0, x] dcost = -0.5 (cosx - 1) = 0.5 (1
      - cosx) cdf(p) = x / 2pi

      theta = acos(1 - 2e1)
      phi = 2pi * e2

      x = sint cosp
      y = sint sinp
      z = cost
    */
    float theta = fm::acos(1 - 2 * sample.x());
    float phi = 2 * fm::pi_f * sample.y();
    if (intersection) {
      float x = fm::sin(theta) * fm::cos(phi);
      float y = fm::sin(theta) * fm::sin(phi);
      float z = fm::cos(theta);

      intersection->position = center + radius * Vector3f{x, y, z};
      intersection->normal = Vector3f{x, y, z};
    }
    if (pdf)
      *pdf = 1 / (4 * fm::pi_f);
  }
  float getArea() const override { return 4 * fm::pi_f * radius * radius; }

public:
  Point3f center;
  float radius;
};