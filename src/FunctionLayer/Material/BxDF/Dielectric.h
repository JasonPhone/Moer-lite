#pragma once
#include "BSDF.h"
#include "Warp.h"

class Dielectric : public BSDF {
 public:
  Dielectric(const Vector3f &_normal, const Vector3f &_tangent,
             const Vector3f &_bitangent, float _ior,
             const Spectrum &_specular_R, const Spectrum &_specular_T)
      : BSDF(_normal, _tangent, _bitangent),
        ior_(_ior),
        inv_ior_(1.f / _ior),
        specular_R_(_specular_R),
        specular_T_(_specular_T) {}

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    return {0.f};
  }

  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const override {
    BSDFSampleResult result;
    auto wo_local = toLocal(wo);
    float eta = wo_local[1] < 0.0 ? ior_ : inv_ior_;

    float cos_theta_T;
    float F = dielectricReflectance(eta, fm::abs(wo_local[1]), &cos_theta_T);
    float reflect_prob = F;
    if (sample[0] <= reflect_prob) {
      result.wi = toWorld(Vector3f{-wo_local[0], wo_local[1], -wo_local[2]});
      result.pdf = reflect_prob;
      result.type = BSDFType::Specular;  // Reflection.
      result.weight = specular_R_ * F / fm::abs(wo_local[1]);
    } else {
      result.wi = toWorld(Vector3f{-eta * wo_local[0],
                                   -std::copysign(cos_theta_T, wo_local[1]),
                                   -eta * wo_local[2]});
      result.pdf = 1 - reflect_prob;
      result.type = BSDFType::Specular;  // Transmission.
      result.weight = specular_T_ * (1 - F) / fm::abs(wo_local[1]);
    }
    return result;
  }

 private:
  float ior_, inv_ior_;
  Spectrum specular_R_, specular_T_;
};