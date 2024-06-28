#pragma once
#include <CoreLayer/Math/Math.h>

inline Vector3f squareToUniformHemisphere(const Vector2f& sample) {
  float y = 1 - 2 * sample[0];
  float r = fm::sqrt(std::max(.0f, 1.f - y * y));
  float phi = 2 * PI * sample[1];
  Vector3f dir{r * fm::sin(phi), std::abs(y), r * fm::cos(phi)};
  return normalize(dir);
}

inline float squareToUniformHemispherePdf(const Vector3f& v) {
  return v[1] >= .0f ? INV_PI * .5f : .0f;
}

inline Vector3f squareToUniformSphere(const Vector2f& sample) {
  float y = 1 - 2 * sample[0];
  float r = fm::sqrt(std::max(.0f, 1.f - y * y));
  float phi = 2 * PI * sample[1];
  Vector3f dir{r * fm::sin(phi), y, r * fm::cos(phi)};
  return normalize(dir);
}

inline float squareToUniformSpherePdf(const Vector3f& v) {
  return INV_PI * .25f;
}

inline Vector3f squareToCosineHemisphere(const Vector2f& sample) {
  float phi = 2 * M_PI * sample[0], theta = std::acos(std::sqrt(sample[1]));
  return Vector3f{std::sin(theta) * std::sin(phi), std::cos(theta),
                  std::sin(theta) * std::cos(phi)};
}

inline float squareToCosineHemispherePdf(const Vector3f& v) {
  return (v[1] > .0f) ? v[1] * INV_PI : .0f;
}

inline float dielectricReflectance(float eta, float cos_theta_I,
                                   float* cos_theta_T) {
  if (cos_theta_I < 0.0) {
    eta = 1.0 / eta;
    cos_theta_I = -cos_theta_I;
  }
  float sin_theta_T_sq = eta * eta * (1.0f - cos_theta_I * cos_theta_I);
  if (sin_theta_T_sq > 1.0) {
    *cos_theta_T = 0.0;
    return 1.0;
  }
  *cos_theta_T = fm::sqrt(std::max(1.0 - sin_theta_T_sq, 0.0));

  float cos_tt = *cos_theta_T;
  float Rs = (eta * cos_theta_I - cos_tt) / (eta * cos_theta_I + cos_tt);
  float Rp = (eta * cos_tt - cos_theta_I) / (eta * cos_tt + cos_theta_I);

  return (Rs * Rs + Rp * Rp) * 0.5;
}