#pragma once
#include <FunctionLayer/Acceleration/Acceleration.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/Mesh.h>

#include "Shape.h"
class TriangleMesh;

class Triangle : public Shape {
 public:
  Triangle() = default;

  Triangle(int _primID, int _vtx0Idx, int _vtx1Idx, int _vtx2Idx,
           const TriangleMesh *_mesh);

  virtual bool rayIntersectShape(Ray &ray, int *primID, float *u,
                                 float *v) const override;

  virtual void fillIntersection(float distance, int primID, float u, float v,
                                Intersection *intersection) const override;

  virtual void uniformSampleOnSurface(Vector2f sample,
                                      Intersection *intersection,
                                      float *pdf) const override {
    // TODO finish this
    return;
  }

 public:
  int primID;
  int vtx0Idx, vtx1Idx, vtx2Idx;
  const TriangleMesh *mesh = nullptr;
};

class TriangleMesh : public Shape {
 public:
  TriangleMesh() = default;

  TriangleMesh(const Json &json);

  //* 当使用embree时，我们使用embree内置的求交函数，故覆盖默认方法
  virtual RTCGeometry getEmbreeGeometry(RTCDevice device) const override;

  virtual bool rayIntersectShape(Ray &ray, int *primID, float *u,
                                 float *v) const override;

  virtual void fillIntersection(float distance, int primID, float u, float v,
                                Intersection *intersection) const override;

  virtual void uniformSampleOnSurface(Vector2f sample,
                                      Intersection *intersection,
                                      float *pdf) const override {
    // TODO finish this
    float u = sample.x();
    int l = 0, r = areaCdf1D.size();
    // Find the interval [a, b): u \in [a, b), which is
    // the index i: u \in [cdf[i], cdf[i+1])
    // Upper bound.

    while (l < r) {
      int mid = (r - l) / 2 + l;
      if (areaCdf1D[mid] <= u)
        l = mid + 1;
      else
        r = mid;
    }
    l--;  // u \in [cdf[l], cdf[l+1])

    u = (u - areaCdf1D[l]) / (areaCdf1D[l + 1] - areaCdf1D[l]);
    fillIntersection(0, l, u, sample.y(), intersection);
    if (pdf) *pdf = (areaCdf1D[l + 1] - areaCdf1D[l]) / getArea(l);
  }

  virtual void initInternalAcceleration() override;

  virtual float getArea() const override;

  float getArea(int faceIdx) const;

  friend class Triangle;

 private:
  std::shared_ptr<MeshData> meshData;
  std::shared_ptr<Acceleration> acceleration;
  std::vector<float> areaCdf1D;
  float area = 0;
};