#pragma once

#include <fstream>

#include "Integrator.h"

namespace RecordKey {
const std::string kLiCalled = "li called";
const std::string kLiBreakByMaxDepth = "li break by max_depth";
const std::string kLiBreakByNoSi = "li break by out of scene";
const std::string kLiContinueByNoBSDF = "li continue by no bsdf";
}  // namespace RecordKey

class SimplePathIntegrator : public Integrator {
 public:
  SimplePathIntegrator() = default;

  SimplePathIntegrator(const Json &json) : Integrator(json) {
    if (json.contains("max_depth") && json["max_depth"].is_number())
      max_depth_ = json["max_depth"];

    if (json.contains("sample_light") && json["sample_light"].is_boolean())
      sample_light_ = json["sample_light"];

    if (json.contains("sample_BSDF") && json["sample_BSDF"].is_boolean())
      sample_BSDF_ = json["sample_BSDF"];

    record.entries[RecordKey::kLiCalled] = 0;
    record.entries[RecordKey::kLiBreakByMaxDepth] = 0;
    record.entries[RecordKey::kLiBreakByNoSi] = 0;
    record.entries[RecordKey::kLiContinueByNoBSDF] = 0;

    log_.open("log.txt", std::ios::out);
  }

  virtual ~SimplePathIntegrator() {
    log_.close();
  }

  virtual Spectrum li(Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

 private:
  int max_depth_ = 1;
  // True: Light is explicitly sampled (uniformly) as end of path.
  //       Infinite lights included.
  // False: Light is sampled only when gets hit.
  bool sample_light_ = true;
  // True: Use BSDF sample to determine next ray direction.
  // False: Next ray goes uniformly.
  bool sample_BSDF_ = true;
  mutable std::ofstream log_;
};
