#pragma once

#include "Integrator.h"

class RandomWalkIntegrator : public Integrator {
 public:
  RandomWalkIntegrator() = default;

  RandomWalkIntegrator(const Json &json) : Integrator(json) {
    if (json["max_depth"].is_number()) {
      max_depth = json["max_depth"];
    } else {
      max_depth = 1;
    }
  }

  virtual ~RandomWalkIntegrator() = default;

  virtual Spectrum li(Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

 private:
  int max_depth;
};
