#include "FunctionLayer/Integrator/RandomWalkIntegrator.h"

#include "CoreLayer/ColorSpace/Spectrum.h"
#include "CoreLayer/Math/Geometry.h"
#include "FunctionLayer/Material/BxDF/BSDF.h"
#include "FunctionLayer/Material/Material.h"
#include "FunctionLayer/Ray/Ray.h"
#include "FunctionLayer/Shape/Intersection.h"

Spectrum RandomWalkIntegrator::li(Ray &ray, const Scene &scene,
                                  std::shared_ptr<Sampler> sampler) const {
  Spectrum spectrum{0.f}, beta{1.f};
  int depth = 0;
  do {
    // Intersection test.
    auto si = scene.rayIntersect(ray);
    if (!si.has_value()) {
      for (auto &light : scene.infiniteLights)
        spectrum += beta * light->evaluateEmission(ray);
      break;
    }
    auto its = si.value();
    computeRayDifferentials(&its, ray);
    // Not evaluating light until hit one.
    if (auto &light = its.shape->light; light)
      spectrum += beta * light->evaluateEmission(its, -ray.direction);
    // Get BSDF at intersection.
    auto material = its.shape->material;
    auto bsdf = material->computeBSDF(its);
    // Sample BSDF randomly.
    auto bs = bsdf->sample(-ray.direction, sampler->next2D());
    // Get leaving ray.
    ray = Ray{its.position, bs.wi};
    beta *= bs.weight * abs(dot(bs.wi, its.normal)) / bs.pdf;

    if (depth++ == max_depth) break;
  } while (1);
  return spectrum;
}

REGISTER_CLASS(RandomWalkIntegrator, "randomWalk")