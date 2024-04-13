#include "AreaLight.h"

#include <ResourceLayer/Factory.h>

#include "FastMath/FastMath.h"
AreaLight::AreaLight(const Json &json) : Light(json) {
  type = LightType::AreaLight;
  shape = Factory::construct_class<Shape>(json["shape"]);

  energy = fetchOptional<Spectrum>(json, "energy", 0.0f);
  power = fetchOptional<Spectrum>(json, "power", 0.0f);
  if (energy.isZero()) {
    // Power(flux) to energy(radiance).
    energy = power / shape->getArea() / fm::pi_f;
  }
}

Spectrum AreaLight::evaluateEmission(const Intersection &intersection,
                                     const Vector3f &wo) const {
  return energy;
}

LightSampleResult AreaLight::sample(const Intersection &shadingPoint,
                                    const Vector2f &sample) const {
  Intersection sampleResult;
  float pdf;
  shape->uniformSampleOnSurface(sample, &sampleResult, &pdf);
  Vector3f shadingPoint2sample = sampleResult.position - shadingPoint.position;

  return {energy,
          normalize(shadingPoint2sample),          // direction
          shadingPoint2sample.length() - EPSILON,  // distance
          sampleResult.normal,                     // normal
          pdf,
          false,
          type};
}

REGISTER_CLASS(AreaLight, "areaLight")