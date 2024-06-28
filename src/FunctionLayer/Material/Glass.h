#pragma once
#include "FunctionLayer/Texture/Texture.h"
#include "Material.h"

class GlassMaterial : public Material {
 public:
  GlassMaterial(const Json &json);

  virtual std::shared_ptr<BSDF> computeBSDF(
      const Intersection &intersection) const override;

 private:
  float ior_;
  std::shared_ptr<Texture<Spectrum>> albedo_R_, albedo_T_;
};
