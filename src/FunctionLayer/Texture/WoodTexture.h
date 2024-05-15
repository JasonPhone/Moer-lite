#pragma once

#include <CoreLayer/ColorSpace/Spectrum.h>

#include "Texture.h"

class WoodTexture : public Texture<Spectrum> {
 public:
  WoodTexture() = delete;
  WoodTexture(const Json& json);
  virtual Spectrum evaluate(const Intersection& intersection) const override;
  virtual Spectrum evaluate(const TextureCoord& texCoord) const override;

 private:
  Vector2f grad(Vector2f sd) const;
  SpectrumRGB c{1.1, .8, .6};
};