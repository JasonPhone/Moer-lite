#pragma once
#include <CoreLayer/Math/Math.h>
//* 使用RGB三通道颜色空间

class SpectrumRGB {
 public:
  SpectrumRGB() : rgb(0.f) {}

  SpectrumRGB(float f) : rgb(f) {}

  SpectrumRGB(float r, float g, float b) : rgb(r, g, b) {}

  SpectrumRGB(Vector3f _rgb) : rgb(_rgb) {}

  SpectrumRGB operator+(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb + rhs.rgb);
  }

  SpectrumRGB &operator+=(const SpectrumRGB &rhs) {
    rgb += rhs.rgb;
    return *this;
  }

  SpectrumRGB operator-(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb - rhs.rgb);
  }

  SpectrumRGB &operator-=(const SpectrumRGB &rhs) {
    rgb -= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator*(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb * rhs.rgb);
  }

  SpectrumRGB &operator*=(const SpectrumRGB &rhs) {
    rgb *= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator*(float f) const { return SpectrumRGB(rgb * f); }

  SpectrumRGB &operator*=(float f) {
    rgb *= f;
    return *this;
  }

  SpectrumRGB operator/(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb / rhs.rgb);
  }

  SpectrumRGB &operator/=(const SpectrumRGB &rhs) {
    rgb /= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator/(float f) const { return SpectrumRGB(rgb / f); }

  SpectrumRGB &operator/=(float f) {
    rgb /= f;
    return *this;
  }

  float operator[](int i) const { return rgb[i]; }

  float &operator[](int i) { return rgb[i]; }

  bool isZero() const { return rgb.isZero(); }

  void debugPrint() const {
    printf("[rgb](");
    for (int i = 0; i < 3; ++i) {
      std::cout << (i == 0 ? '\0' : ',') << rgb[i];
    }
    printf(")%c", '\n');
    fflush(stdout);
  }

  float Y() const {
    /**
     * @brief Adobe RGB(1998) to XYZ in D50,
     *        see http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
     */
    return rgb[0] * 0.3111242 + rgb[1] * 0.6256560 + rgb[2] * 0.0632197;
  }

 private:
  Vector3f rgb;
};

inline SpectrumRGB operator*(float f, const SpectrumRGB &spectrum) {
  return spectrum * f;
}

inline Vector3f toVec3(const SpectrumRGB &spectrum) {
  return Vector3f{spectrum[0], spectrum[1], spectrum[2]};
}

inline SpectrumRGB toSpectrum(const Vector3f &vec) {
  return SpectrumRGB(vec[0], vec[1], vec[2]);
}