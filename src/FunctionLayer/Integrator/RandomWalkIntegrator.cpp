#include "FunctionLayer/Integrator/RandomWalkIntegrator.h"
#include "CoreLayer/ColorSpace/Spectrum.h"
#include "CoreLayer/Math/Geometry.h"
#include "FunctionLayer/Material/BxDF/BSDF.h"
#include "FunctionLayer/Material/Material.h"
#include "FunctionLayer/Ray/Ray.h"
#include "FunctionLayer/Shape/Intersection.h"

Spectrum RandomWalkIntegrator::li(Ray &ray, const Scene &scene,
                                  std::shared_ptr<Sampler> sampler) const {
  return liRandomWalk(ray, scene, sampler, 0);
}
Spectrum RandomWalkIntegrator::liRandomWalk(Ray &ray, const Scene &scene,
                                            std::shared_ptr<Sampler> sampler,
                                            int depth) const {
  Spectrum Linf{0.f};
  // Intersection test.
  auto intersection_opt = scene.rayIntersect(ray);
  if (!intersection_opt.has_value()) {
    for (auto &light : scene.infiniteLights)
      Linf += light->evaluateEmission(ray);
    return Linf;
  }
  auto intersection = intersection_opt.value();
  computeRayDifferentials(&intersection, ray);
  // Le.
  Spectrum Le{0.f};
  if (auto &light = intersection.shape->light; light)
    Le += light->evaluateEmission(intersection, -ray.direction);
  // If should terminate.
  if (depth == max_depth)
    return Le;
  // Get BSDF at intersection.
  auto material = intersection.shape->material;
  auto bsdf = material->computeBSDF(intersection);
  // Sample BSDF.
  auto bsdf_sample = bsdf->sample(-ray.direction, sampler->next2D());
  Spectrum f{0.f};
  // Get leaving ray.
  auto next_ray = Ray{intersection.position, bsdf_sample.wi};
  // Recurse for next step.
  f = liRandomWalk(next_ray, scene, sampler, depth + 1);
  f *= bsdf_sample.weight / bsdf_sample.pdf;

  return Le + f;
}

REGISTER_CLASS(RandomWalkIntegrator, "randomWalk")