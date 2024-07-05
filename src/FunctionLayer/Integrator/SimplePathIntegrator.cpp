#include "FunctionLayer/Integrator/SimplePathIntegrator.h"

#include "CoreLayer/ColorSpace/Spectrum.h"
#include "CoreLayer/Math/Geometry.h"
#include "FunctionLayer/Material/BxDF/BSDF.h"
#include "FunctionLayer/Material/BxDF/Warp.h"
#include "FunctionLayer/Material/Material.h"
#include "FunctionLayer/Ray/Ray.h"
#include "FunctionLayer/Shape/Intersection.h"

Spectrum SimplePathIntegrator::li(Ray &ray, const Scene &scene,
                                  std::shared_ptr<Sampler> sampler) const {
  record.entries[RecordKey::kLiCalled]++;
  Spectrum L{0.f}, beta{1.f};
  int depth = 0;
  // True if current ray is from specular BSDF.
  bool from_specular = true;
  do {
    // Test along ray.
    auto si = scene.rayIntersect(ray);
    if (!si.has_value()) {
      if (!sample_light_ || from_specular)
        // Lights are not sampled until hit,
        // so infinite lights should be added, or they'll never get sampled.
        // Lights are also not sampled if last BSDF is specular,
        // under that condition the direction is fixed.
        for (auto &light : scene.infiniteLights)
          L += beta * light->evaluateEmission(ray);
      record.entries[RecordKey::kLiBreakByNoSi]++;
      break;
    }

    auto its = si.value();
    Vector3f wo = -ray.direction;

    if (!sample_light_ || from_specular)
      // from_specular also enables direct hit on light.
      if (auto light = its.shape->light; light)
        L += beta * light->evaluateEmission(its, wo);

    if (depth++ == max_depth_) {
      record.entries[RecordKey::kLiBreakByMaxDepth]++;
      break;
    }
    // Get BSDF.
    auto bsdf = its.shape->material->computeBSDF(its);
    if (!bsdf) {
      // Skip to next intersection in medium rendering.
      ray = Ray(its.position + EPSILON * ray.direction, ray.direction);
      from_specular = true;
      record.entries[RecordKey::kLiContinueByNoBSDF]++;
      continue;
    }
    // Sample lights as p_{end}.
    if (sample_light_) {
      // Infinite lights.
      for (auto light : scene.infiniteLights) {
        auto light_sample = light->sample(its, sampler->next2D());
        Ray shadow_ray(its.position, light_sample.direction, 1e-4f, FLT_MAX);
        // Light connected.
        if (auto occlude = scene.rayIntersect(shadow_ray);
            !occlude.has_value()) {
          // Spectrum f = bsdf->f(wo, shadow_ray.direction);
          Spectrum f = bsdf->f(wo, shadow_ray.direction) *
                       abs(dot(shadow_ray.direction, its.normal));
          float pdf = convertPDF(light_sample, its);
          L += beta * light_sample.energy * f / pdf;
        }
      }

      float pdf_light = .0f;
      // Common lights.
      auto light = scene.sampleLight(sampler->next1D(), &pdf_light);
      if (light && pdf_light != .0f) {
        auto light_sample = light->sample(its, sampler->next2D());
        Ray shadow_ray(its.position, light_sample.direction, 1e-4f,
                       light_sample.distance);
        if (auto occlude = scene.rayIntersect(shadow_ray);
            !occlude.has_value()) {
          // Spectrum f = bsdf->f(wo, shadow_ray.direction);
          Spectrum f = bsdf->f(wo, shadow_ray.direction) *
                       abs(dot(shadow_ray.direction, its.normal));
          light_sample.pdf *= pdf_light;
          float pdf = convertPDF(light_sample, its);
          L += beta * light_sample.energy * f / pdf;
        }
      }
    }
    // Current at p_i, query for p_{i+1}. Sample BSDF or sphere surface.
    // Construct path, accumulate throughput.
    Vector3f wi;
    if (sample_BSDF_) {
      auto bs = bsdf->sample(wo, sampler->next2D());
      // BSDF out direction.
      // The cos term is from geometry.
      beta *= bs.weight * abs(dot(bs.wi, its.normal)) / bs.pdf;
      // beta *= bs.weight * abs(dot(bs.wi, its.normal));
      from_specular = bs.type == BSDFType::Specular;
      wi = bs.wi;
    } else {
      // Uniform out direction.
      // TODO Aggregated BSDF types for better direction correction.
      wi = squareToUniformSphere(sampler->next2D());
      float pdf = squareToUniformSpherePdf(wi);
      beta *= bsdf->f(wo, wi) * abs(dot(wi, its.normal)) / pdf;
      // beta *= bsdf->f(wo, wi) * abs(dot(wi, its.normal));
      from_specular = false;  // Note the correlation.
    }
    ray = Ray{its.position + EPSILON * wi, wi};
  } while (!beta.isZero());  // Early stop.
  return L;
}

REGISTER_CLASS(SimplePathIntegrator, "simplePath")