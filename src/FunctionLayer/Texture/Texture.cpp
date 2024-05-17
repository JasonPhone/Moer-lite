#include "Texture.h"

TextureCoord UVMapping::map(const Intersection &intersection) const {
  Vector2f coord = intersection.texCoord;
  /// @brief Sphere texture mapping.
  // Vector3f np = normalize(intersection.position - Point3f{0, 0, 0});
  // float phi = fm::acos(np[2]), theta = fm::atan(np[1] / np[0]);
  // float u = theta * 0.5 * INV_PI, v = phi * INV_PI;
  /// @brief Cylinder texture mapping.
  // Vector3f np = normalize(intersection.position - Point3f{0, 0, 0});
  // float phi = fm::acos(np[2]), theta = fm::atan(np[1] / np[0]);
  // float u = theta * 0.5 * INV_PI, v = np[2] >= 0 ? np[2] : np[2] + 1;
  /// @brief Flat texture mapping.
  // Vector3f v0{1, 0, 0}, v1{0, 1, 0}, p = intersection.position - Point3f{0, 0, 0};
  // float u = dot(v0, p), v = dot(v1, p);

  // coord[0] = u - int(u), coord[1] = v - int(v);
  // coord[0] += coord[0] < 0;
  // coord[1] += coord[1] < 0;
  return TextureCoord{coord,
                      Vector2f{intersection.dudx, intersection.dvdx},
                      Vector2f{intersection.dudy, intersection.dvdy}};
}