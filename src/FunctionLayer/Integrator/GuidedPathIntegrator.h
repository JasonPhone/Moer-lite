#pragma once

#include <fstream>

#include "Integrator.h"

class STree;

struct RadianceRecord {
  RadianceRecord(const Vector3f &d, float r, float p)
      : dir(d), radiance(r), pdf(p) {}
  Vector3f dir;
  float radiance;
  // Pdf w.r.t to current vertex solid angle domain.
  float pdf;
};
// The closer one to camera of two path vertices.
struct PathVertex {
  PathVertex(const Point3f &p, float _pdf) : pos{p}, pdf(_pdf) {}
  Point3f pos;
  float pdf;
  std::vector<RadianceRecord> records;
};

class GuidedPathIntegrator : public Integrator {
 public:
  GuidedPathIntegrator() = delete;

  GuidedPathIntegrator(const Json &json);

  virtual ~GuidedPathIntegrator();

  virtual void adaptScene(const Scene &scene) override;
  virtual void preSampling() override;
  virtual void postSampling(int spp) override;

  virtual Spectrum li(Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

 private:
  int m_max_depth = 1;
  // True: Light is explicitly sampled (uniformly) as end of path.
  //       Infinite lights included.
  // False: Light is sampled only when gets hit.
  bool m_sample_light = true;
  mutable std::ofstream m_log;
  int m_n_iteration = 0;
  STree *m_sampling_STree = nullptr, *m_building_STree = nullptr;
};
