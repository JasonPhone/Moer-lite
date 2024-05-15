#include "WoodTexture.h"

#include <ResourceLayer/Factory.h>

WoodTexture::WoodTexture(const Json& json) : Texture<Spectrum>() {}

Spectrum WoodTexture::evaluate(const Intersection& intersection) const {
  return evaluate(TextureCoord{intersection.texCoord});
}
Spectrum WoodTexture::evaluate(const TextureCoord& texCoord) const {
  Vector2f p = texCoord.coord;

  p.a[0] *= 6, p.a[1] *= 6;
  Vector2f l[4], v[4], d[4];
  float f[4];
  for (int i = 0; i < 4; i++) {
    float a = int(p.a[0] + (i & 1));
    float b = int(p.a[1] + (i & 2));

    l[i] = Vector2f{a, b};
    v[i] = grad(l[i]);
    d[i] = p - v[i];
    f[i] = dot(d[i], v[i]);
  }

  float rx = (p.a[0] - l[0].a[0]) / (l[1].a[0] - l[0].a[0]);
  float ry = (p.a[1] - l[0].a[1]) / (l[2].a[1] - l[0].a[1]);

  float x_low = (1 - rx) * f[0] + rx * f[1];
  float x_high = (1 - rx) * f[2] + rx * f[3];
  float sample = (1 - ry) * x_low + ry * x_high;

  return (sample * 5 - int(sample * 5)) * c;
}

Vector2f WoodTexture::grad(Vector2f sd) const {
  uint32_t x = sd.a[0], y = sd.a[1];
  uint32_t s = (x << 16) | y;
  std::srand(s * s);
  Vector2f g{1.f * std::rand(), 1.f * std::rand()};
  g /= g.len();
  return g;
}

REGISTER_CLASS(WoodTexture, "woodTex")