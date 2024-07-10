#include "Glass.h"

#include "./BxDF/Dielectric.h"
#include "FunctionLayer/Texture/ConstantTexture.h"

GlassMaterial::GlassMaterial(const Json& json) : Material(json) {
  // ior, R, T.
  ior_ = fetchOptional(json, "ior", 1.52);  // Common glass.
  auto s = fetchRequired<Spectrum>(json, "albedoR");
  albedo_R_ = std::make_shared<ConstantTexture<Spectrum>>(s);
  s = fetchRequired<Spectrum>(json, "albedoT");
  albedo_T_ = std::make_shared<ConstantTexture<Spectrum>>(s);
}

std::shared_ptr<BSDF> GlassMaterial::computeBSDF(
    const Intersection& intersection) const {
  Vector3f normal, tangent, bitangent;
  computeShadingGeometry(intersection, &normal, &tangent, &bitangent);
  Spectrum spec_R = albedo_R_->evaluate(intersection);
  Spectrum spec_T = albedo_T_->evaluate(intersection);
  return std::make_shared<Dielectric>(normal, tangent, bitangent, ior_, spec_R,
                                      spec_T);
}

REGISTER_CLASS(GlassMaterial, "glass")