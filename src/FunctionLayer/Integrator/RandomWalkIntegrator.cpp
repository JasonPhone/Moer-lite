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
    Spectrum Linf{0.f};
    // Intersection test.
    auto intersection_opt = scene.rayIntersect(ray);
    if (!intersection_opt.has_value()) {
      for (auto &light : scene.infiniteLights)
        spectrum += beta * light->evaluateEmission(ray);
      break;
    }
    auto intersection = intersection_opt.value();
    computeRayDifferentials(&intersection, ray);
    // Le.
    if (auto &light = intersection.shape->light; light)
      spectrum += beta * light->evaluateEmission(intersection, -ray.direction);
    // If should terminate.
    if (depth == max_depth)
      break;
    // Get BSDF at intersection.
    auto material = intersection.shape->material;
    auto bsdf = material->computeBSDF(intersection);
    // Sample BSDF.
    auto bsdf_sample = bsdf->sample(-ray.direction, sampler->next2D());
    // Get leaving ray.
    ray = Ray{intersection.position, bsdf_sample.wi};
    beta *= bsdf_sample.weight;
    // Next step.
  } while (1);
  return spectrum;
}
Spectrum RandomWalkIntegrator::liRandomWalk(Ray &ray, const Scene &scene,
                                            std::shared_ptr<Sampler> sampler,
                                            int depth) const {}

REGISTER_CLASS(RandomWalkIntegrator, "randomWalk")