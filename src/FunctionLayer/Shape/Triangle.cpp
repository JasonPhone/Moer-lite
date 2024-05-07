#include "Triangle.h"

#include <FunctionLayer/Acceleration/Linear.h>
//--- Triangle ---
Triangle::Triangle(int _primID, int _vtx0Idx, int _vtx1Idx, int _vtx2Idx,
                   const TriangleMesh *_mesh)
    : primID(_primID),
      vtx0Idx(_vtx0Idx),
      vtx1Idx(_vtx1Idx),
      vtx2Idx(_vtx2Idx),
      mesh(_mesh) {
  Point3f vtx0 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx0Idx]),
          vtx1 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx1Idx]),
          vtx2 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx2Idx]);
  boundingBox.Expand(vtx0);
  boundingBox.Expand(vtx1);
  boundingBox.Expand(vtx2);
  this->geometryID = mesh->geometryID;
}

bool Triangle::rayIntersectShape(Ray &ray, int *primID, float *u,
                                 float *v) const {
  Point3f origin = ray.origin;
  Vector3f direction = ray.direction;
  Point3f vtx0 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx0Idx]),
          vtx1 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx1Idx]),
          vtx2 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx2Idx]);

  Vector3f edge0 = vtx1 - vtx0, edge1 = vtx2 - vtx0;

  Vector3f paralNormal = normalize(cross(edge0, edge1));
  float d = -dot(paralNormal, Vector3f{vtx0[0], vtx0[1], vtx0[2]});
  float a = dot(paralNormal, Vector3f{origin[0], origin[1], origin[2]}) + d;
  float b = dot(paralNormal, direction);
  if (b == .0f) return false;  // miss
  float t = -a / b;

  if (t < ray.tNear || t > ray.tFar) return false;

  Point3f hitpoint = origin + t * direction;
  // hitpoint = vtx0 + u * e0 + v * e1, 0 <= u, v <= 1
  Vector3f v1 = cross(hitpoint - vtx0, edge1), v2 = cross(edge0, edge1);
  float u_ = v1.length() / v2.length();
  if (dot(v1, v2) < 0) u_ *= -1;

  v1 = cross(hitpoint - vtx0, edge0);
  v2 = cross(edge1, edge0);
  float v_ = v1.length() / v2.length();
  if (dot(v1, v2) < 0) v_ *= -1;

  if (u_ >= .0f && v_ >= .0f && (u_ + v_ <= 1.f)) {
    ray.tFar = t;
    *primID = this->primID;
    *u = u_;
    *v = v_;
    return true;
  }

  return false;
}

void Triangle::fillIntersection(float distance, int primID, float u, float v,
                                Intersection *intersection) const {
  // 该函数实际上不会被调用
  return;
}

// float Triangle::getArea() const {
//   Point3f v0 = mesh->meshData->vertexBuffer[vtx0Idx];
//   Point3f v1 = mesh->meshData->vertexBuffer[vtx1Idx];
//   Point3f v2 = mesh->meshData->vertexBuffer[vtx2Idx];
//   Vector3f v01 = v1 - v0;
//   Vector3f v02 = v2 - v0;
//   return cross(v01, v02).length() * 0.5;
// }
//--- TriangleMesh ---
TriangleMesh::TriangleMesh(const Json &json) : Shape(json) {
  const auto &filepath = fetchRequired<std::string>(json, "file");
  meshData = MeshData::loadFromFile(filepath);
  // For triangle sampling based on area portion.
  areaCdf1D.push_back(0);
  for (int i = 0; i < meshData->faceCount; i++) {
    float a = getArea(i);
    areaCdf1D.push_back(a);
    area += a;
  }
  for (int i = 1; i < areaCdf1D.size(); i++) {
    areaCdf1D[i] /= area;
    areaCdf1D[i] += areaCdf1D[i - 1];
  }
}

RTCGeometry TriangleMesh::getEmbreeGeometry(RTCDevice device) const {
  RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

  float *vertexBuffer = (float *)rtcSetNewGeometryBuffer(
      geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float),
      meshData->vertexCount);
  for (int i = 0; i < meshData->vertexCount; ++i) {
    Point3f vertex = transform.toWorld(meshData->vertexBuffer[i]);
    vertexBuffer[3 * i] = vertex[0];
    vertexBuffer[3 * i + 1] = vertex[1];
    vertexBuffer[3 * i + 2] = vertex[2];
  }

  unsigned *indexBuffer = (unsigned *)rtcSetNewGeometryBuffer(
      geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
      3 * sizeof(unsigned), meshData->faceCount);
  for (int i = 0; i < meshData->faceCount; ++i) {
    indexBuffer[i * 3] = meshData->faceBuffer[i][0].vertexIndex;
    indexBuffer[i * 3 + 1] = meshData->faceBuffer[i][1].vertexIndex;
    indexBuffer[i * 3 + 2] = meshData->faceBuffer[i][2].vertexIndex;
  }
  rtcCommitGeometry(geometry);
  return geometry;
}

bool TriangleMesh::rayIntersectShape(Ray &ray, int *primID, float *u,
                                     float *v) const {
  //* 当使用embree加速时，该方法不会被调用
  int geomID = -1;
  return acceleration->rayIntersect(ray, &geomID, primID, u, v);
}

void TriangleMesh::fillIntersection(float distance, int primID, float u,
                                    float v, Intersection *intersection) const {
  intersection->distance = distance;
  intersection->shape = this;
  //* 在三角形内部用插值计算交点、法线以及纹理坐标
  auto faceInfo = meshData->faceBuffer[primID];
  float w = 1.f - u - v;

  //* 计算交点
  Point3f pw = transform.toWorld(
              meshData->vertexBuffer[faceInfo[0].vertexIndex]),
          pu = transform.toWorld(
              meshData->vertexBuffer[faceInfo[1].vertexIndex]),
          pv = transform.toWorld(
              meshData->vertexBuffer[faceInfo[2].vertexIndex]);
  intersection->position = Point3f{w * pw[0] + u * pu[0] + v * pv[0],
                                   w * pw[1] + u * pu[1] + v * pv[1],
                                   w * pw[2] + u * pu[2] + v * pv[2]};
  //* 计算法线
  if (meshData->normalBuffer.size() != 0) {
    Vector3f nw = transform.toWorld(
                 meshData->normalBuffer[faceInfo[0].normalIndex]),
             nu = transform.toWorld(
                 meshData->normalBuffer[faceInfo[1].normalIndex]),
             nv = transform.toWorld(
                 meshData->normalBuffer[faceInfo[2].normalIndex]);
    intersection->normal = normalize(w * nw + u * nu + v * nv);
  } else {
    intersection->normal = normalize(cross(pu - pw, pv - pw));
  }

  //* 计算纹理坐标
  if (meshData->texcodBuffer.size() != 0) {
    Vector2f tw = meshData->texcodBuffer[faceInfo[0].texcodIndex],
             tu = meshData->texcodBuffer[faceInfo[1].texcodIndex],
             tv = meshData->texcodBuffer[faceInfo[2].texcodIndex];
    intersection->texCoord = w * tw + u * tu + v * tv;
  } else {
    intersection->texCoord = Vector2f{.0f, .0f};
  }

  // 计算交点的切线和副切线
  Vector3f tangent{1.f, 0.f, .0f};
  Vector3f bitangent;
  if (std::abs(dot(tangent, intersection->normal)) > .9f) {
    tangent = Vector3f(.0f, 1.f, .0f);
  }
  bitangent = normalize(cross(tangent, intersection->normal));
  tangent = normalize(cross(intersection->normal, bitangent));
  intersection->tangent = tangent;
  intersection->bitangent = bitangent;

  // dpdu and dpdv.
  /*
    p1 = p0 + Du01 dpdu + Dv01 dpdv
    p2 = p0 + Du02 dpdu + Dv02 dpdv
    dp01 = Du01 dpdu + Dv01 dpdv
    dp02 = Du02 dpdu + Dv02 dpdv

    dp01 = Du01 Dv01 dot dpdu
    dp02   Du02 Dv02     dpdv
  */
  Vector3f dpdu, dpdv;
  Point3f p0 = transform.toWorld(
              meshData->vertexBuffer[faceInfo[0].vertexIndex]),
          p1 = transform.toWorld(
              meshData->vertexBuffer[faceInfo[1].vertexIndex]),
          p2 = transform.toWorld(
              meshData->vertexBuffer[faceInfo[2].vertexIndex]);
  Vector2f t0, t1, t2;

  if (meshData->texcodBuffer.size() != 0) {
    t0 = meshData->texcodBuffer[faceInfo[0].texcodIndex];
    t1 = meshData->texcodBuffer[faceInfo[1].texcodIndex];
    t2 = meshData->texcodBuffer[faceInfo[2].texcodIndex];
  } else {
    t0 = Vector2f{0, 0};
    t1 = Vector2f{0, 1};
    t2 = Vector2f{1, 0};
  }
  Vector3f dp01 = p1 - p0, dp02 = p2 - p0;
  Vector2f dt01 = t1 - t0, dt02 = t2 - t0;
  /*
    dp01^T = dt01^T dot dpdu^T
    dp02^T   dt02^T     dpdv^T
  */
  float det = dt01[0] * dt02[1] - dt01[1] * dt02[0];
  if (std::abs(det) >= 1e-8) {
    float invdet = 1.0 / det;
    dpdu = invdet * (dt02[1] * dp01 - dt01[1] * dp02);
    dpdv = invdet * (dt01[0] * dp02 - dt02[0] * dp01);
  }

  intersection->dpdu = dpdu;
  intersection->dpdv = dpdv;
}

void TriangleMesh::initInternalAcceleration() {
  acceleration = Acceleration::createAcceleration();
  int primCount = meshData->faceCount;
  for (int primID = 0; primID < primCount; ++primID) {
    int vtx0Idx = meshData->faceBuffer[primID][0].vertexIndex,
        vtx1Idx = meshData->faceBuffer[primID][1].vertexIndex,
        vtx2Idx = meshData->faceBuffer[primID][2].vertexIndex;
    std::shared_ptr<Triangle> triangle =
        std::make_shared<Triangle>(primID, vtx0Idx, vtx1Idx, vtx2Idx, this);
    acceleration->attachShape(triangle);
  }
  acceleration->build();
  // TriangleMesh的包围盒就是其内部加速结构的包围盒
  boundingBox = acceleration->boundingBox;
}

float TriangleMesh::getArea() const { return area; }

float TriangleMesh::getArea(int faceIdx) const {
  auto faceData = meshData->faceBuffer[faceIdx];

  Point3f p0 = transform.toWorld(
              meshData->vertexBuffer[faceData[0].vertexIndex]),
          p1 = transform.toWorld(
              meshData->vertexBuffer[faceData[1].vertexIndex]),
          p2 = transform.toWorld(
              meshData->vertexBuffer[faceData[2].vertexIndex]);

  Vector3f v01 = p1 - p0;
  Vector3f v02 = p2 - p0;

  return cross(v01, v02).length() * 0.5;
}

REGISTER_CLASS(TriangleMesh, "triangle")