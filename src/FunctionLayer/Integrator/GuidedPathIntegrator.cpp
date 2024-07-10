#include "FunctionLayer/Integrator/GuidedPathIntegrator.h"

#include <stack>

#include "CoreLayer/ColorSpace/Spectrum.h"
#include "CoreLayer/Math/Geometry.h"
#include "FunctionLayer/Material/BxDF/BSDF.h"
#include "FunctionLayer/Material/BxDF/Warp.h"
#include "FunctionLayer/Material/Material.h"
#include "FunctionLayer/Ray/Ray.h"
#include "FunctionLayer/Shape/Intersection.h"

// All in world space.

struct DTreeNode {
  bool is_leaf = true;
  float children_flux[4] = {0, 0, 0, 0};
  int children_idx[4] = {-1, -1, -1, -1};
};

struct SDTreeSampleResult {
  Vector3f wi;
  float pdf;
};

/// @brief {theta, phi} in [0, 1)^2 to {x, y, z} on sphere surface.
Vector3f canonicalToDir(const Vector2f &p) {
  float theta = p[0], phi = p[1];
  theta = clamp(theta, 0.f + EPSILON, 1.f - EPSILON);
  phi = clamp(phi, 0.f + EPSILON, 1.f - EPSILON);
  theta *= PI;
  phi = (phi * 2 * PI) - PI;
  float cos_t = fm::cos(theta), sin_t = fm::sin(theta);
  float cos_p = fm::cos(phi), sin_p = fm::sin(phi);
  return {sin_p * sin_t, cos_t, cos_p * sin_t};
}

/// @brief {x, y, z} on sphere surface to {theta, phi} in [0, 1)^2.
Vector2f dirToCanonical(const Vector3f &d) {
  float sinp_sint = d[0], cos_t = d[1], cosp_sint = d[2];
  float theta = fm::acos(cos_t);                // [0, pi].
  float phi = fm::atan2(sinp_sint, cosp_sint);  // [-pi, +pi].

  theta *= INV_PI;
  phi = (phi + PI) * 0.5f * INV_PI;
  theta = clamp(theta, 0.f + EPSILON, 1.f - EPSILON);
  phi = clamp(phi, 0.f + EPSILON, 1.f - EPSILON);
  return {theta, phi};
}

class DTree {
  /**
   * Spherical Coord.
   * phi [0, 2pi)
   * |________
   * |10 |11 |
   * |___|___|
   * |00 |01 |
   * |___|___|__ theta [0, pi)
   *
   */
 public:
  DTree() {
    m_nodes.clear();
    m_nodes.emplace_back();
    m_nodes[0].is_leaf = true;
    std::fill(m_nodes[0].children_flux, m_nodes[0].children_flux + 4, 0);
    std::fill(m_nodes[0].children_idx, m_nodes[0].children_idx + 4, -1);
  }
  ~DTree() {}
  void clear() {
    m_nodes.clear();
    m_nodes.emplace_back();
    m_nodes[0].is_leaf = true;
    std::fill(m_nodes[0].children_flux, m_nodes[0].children_flux + 4, 0);
    std::fill(m_nodes[0].children_idx, m_nodes[0].children_idx + 4, -1);
  }
  void sample(const Vector2f &sample, SDTreeSampleResult &sds) const {
    // sds.wi = squareToCosineHemisphere(sample);
    // sds.pdf = squareToCosineHemispherePdf(sds.wi);
    // return;
    // NOTE Sample randomness may lost if tree goes too deep.
    int cur_node_idx = 0;
    float e_theta = sample[0], e_phi = sample[1];
    float theta = 0, phi = 0;  // Normalized [0, 1).
    float range_scale = 1.f;   // Angle extent a node covers, shared.
    // Importance sampling, pdf should be flux(this node) / flux(entire tree).
    float pdf = 0.25 * PI;
    while (!m_nodes[cur_node_idx].is_leaf) {
      const auto &cur_node = m_nodes[cur_node_idx];
      float flux00 = cur_node.children_flux[0];
      float flux01 = cur_node.children_flux[1];
      float flux10 = cur_node.children_flux[2];
      float flux11 = cur_node.children_flux[3];
      float flux_sum = flux00 + flux01 + flux10 + flux11;
      // Theta first, then phi.
      float sample_thresh = (flux00 + flux10) / flux_sum;
      int child_idx = 0;
      if (e_theta < sample_thresh) {
        // Sample 00+10.
        e_theta /= sample_thresh;                    // Restore randomness.
        child_idx += 0b00;                           // Select child.
        theta += 0 * range_scale;                    // Update angle.
        sample_thresh = flux00 / (flux00 + flux10);  // Update threshold.
      } else {
        // Sample 01+11.
        e_theta = (e_theta - sample_thresh) /
                  (1.f - sample_thresh);             // Restore randomness.
        child_idx += 0b01;                           // Select child.
        theta += 0.5f * range_scale;                 // Update angle.
        sample_thresh = flux01 / (flux01 + flux11);  // Update threshold.
      }
      if (e_phi < sample_thresh) {
        // Sample 00 or 01.
        e_phi /= sample_thresh;  // Restore randomness.
        child_idx += 0b00;       // Select child.
        phi += 0 * range_scale;  // Update angle.
      } else {
        // Sample 10 or 11.
        e_phi = (e_phi - sample_thresh) /
                (1.f - sample_thresh);  // Restore randomness.
        child_idx += 0b10;              // Select child.
        phi += 0.5f * range_scale;      // Update angle.
      }

      pdf *= 4 * cur_node.children_flux[child_idx] / flux_sum;

      range_scale *= 0.5f;
      cur_node_idx = cur_node.children_idx[child_idx];
    }

    // SDTreeSampleResult res;
    // Random select in current range.
    theta += range_scale * e_theta;
    phi += range_scale * e_phi;
    sds.wi = canonicalToDir(Vector2f{theta, phi});

    // sds.pdf = pdf;
    // pdf = pdf * range_scale * range_scale;
    sds.pdf = pdf * abs(sin(theta));
  }

  void update(const Vector3f &wi, float irradiance) {
    Vector3f w = normalize(wi);
    auto ca = dirToCanonical(w);
    auto di = canonicalToDir(ca);
    for (int i = 0; i < 3; i++) {
      float diff = abs(wi[i] - di[i]);
      if (diff > EPSILON * 10) {
        w.debugPrint();
        ca.debug_print();
        di.debugPrint();
        assert(false);
      }
    }
    float theta = ca[0], phi = ca[1];
    float scale = 1;
    int cur_node_idx = 0;
    while (!m_nodes[cur_node_idx].is_leaf) {
      auto &cur_node = m_nodes[cur_node_idx];
      int child_index = 0;
      if (theta < 0.5f) {
        child_index += 0b00;
        theta *= 2;
      } else {
        child_index += 0b01;
        theta = (theta - 0.5f) * 2;
      }
      if (phi < 0.5f) {
        child_index += 0b00;
        phi *= 2;
      } else {
        child_index += 0b10;
        phi = (phi - 0.5f) * 2;
      }
      cur_node.children_flux[child_index] += irradiance;
      cur_node_idx = cur_node.children_idx[child_index];
    }
    for (int i = 0; i < 4; i++)
      m_nodes[cur_node_idx].children_flux[i] += irradiance;
  }

  void refine() {
    float total_flux = nodeFlux(m_nodes[0]);
    if (total_flux == 0) return;
    std::stack<int> nodes;
    nodes.push(0);
    while (!nodes.empty()) {
      int at_idx = nodes.top();
      nodes.pop();
      if (at_idx < 0) continue;
      DTreeNode &at = m_nodes[at_idx];
      float node_flux = nodeFlux(at);
      if (node_flux / total_flux <= m_flux_balance)
        // prune(at);
        ;
      else if (at.is_leaf)
        subdivide(at);
      for (int i = 0; i < 4; i++) nodes.push(at.children_idx[i]);
    }
    // std::cout << "DTree refine end, #nodes " << m_nodes.size()
    //           << ", total flux " << nodeFlux(m_nodes[0]) << std::endl;
  }

  void debugPrint() {
    // for (int i = 0; i < m_nodes.size(); i++) {
    //   auto &node = m_nodes[i];
    //   std::cout << "node\t" << i << ":\n";
    //   std::cout << "\tflux\t" << node.children_flux[0] << ", "
    //             << node.children_flux[1] << ", " << node.children_flux[2]
    //             << ", " << node.children_flux[3] << "\n";
    //   std::cout << "\tchild\t" << node.children_idx[0] << ", "
    //             << node.children_idx[1] << ", " << node.children_idx[2] << ",
    //             "
    //             << node.children_idx[3] << "\n";
    // }
  }

  int getNumSamples() const { return m_n_samples; }
  void setNumSamples(int val) { m_n_samples = val; }

 private:
  void subdivide(DTreeNode &node) {
    float children_flux = nodeFlux(node) * 0.25f * 0.25f;
    for (int i = 0; i < 4; i++) {
      m_nodes.emplace_back();
      auto &new_node = m_nodes.back();
      new_node.is_leaf = true;
      std::fill(new_node.children_flux, new_node.children_flux + 4,
                children_flux);
      std::fill(new_node.children_idx, new_node.children_idx + 4, -1);
      node.children_idx[i] = m_nodes.size() - 1;
    }
    node.is_leaf = false;
  }
  // void prune(DTreeNode *node) {
  //   if (node->is_leaf) return;
  //   deleteChildren(node);
  // }
  // void deleteChildren(DTreeNode *node) {
  //   if (node == nullptr) return;
  //   for (int i = 0; i < 4; i++) {
  //     deleteChildren(node->children[i]);
  //     delete node->children[i];
  //     node->children[i] = nullptr;
  //   }
  //   node->is_leaf = true;
  // }
  float nodeFlux(const DTreeNode &node) const {
    return node.children_flux[0] + node.children_flux[1] +
           node.children_flux[2] + node.children_flux[3];
  }
  int m_n_samples = 0;
  float m_flux_balance = 0.01f;
  std::vector<DTreeNode> m_nodes;
};

struct STreeNode {
  DTree d_tree;  // Number of samples are stored in corresponding DTree.
  int axis_idx = 0;
  bool is_leaf = true;
  int children_idx[2] = {-1, -1};

  float bvh_extend = 0;
};

class STree {
 public:
  STree(const AABB &aabb) : m_aabb{aabb} {
    m_nodes.clear();
    m_nodes.emplace_back();  // Empty root.
    m_nodes[0].is_leaf = true;

    // Extend AABB to a full cube.
    Vector3f sz = m_aabb.pMax - m_aabb.pMin;
    float max_extend = std::max(sz[0], std::max(sz[1], sz[2]));
    m_aabb.pMax = m_aabb.pMin + Vector3f(max_extend);

    m_nodes[0].bvh_extend = max_extend;
  }
  void clear() {
    for (auto &node : m_nodes) {
      if (node.is_leaf) node.d_tree.clear();
    }
  }
  void refine(int last_itr_spp) {
    assert(m_nodes.size() < std::numeric_limits<size_t>::max() - 1);
    int threshold = last_itr_spp;

    struct StkNode {
      int index, depth;
    };
    std::stack<StkNode> nodes;
    nodes.push({0, 1});
    while (!nodes.empty()) {
      auto node = nodes.top();
      nodes.pop();
      if (node.index < 0 || node.index >= m_nodes.size()) continue;
      if (m_nodes[node.index].is_leaf) {
        if (m_nodes[node.index].d_tree.getNumSamples() > threshold) {
          subdivide(node.index);
        }
      } else {
        const auto &p = m_nodes[node.index];
        nodes.push({p.children_idx[0], node.depth + 1});
        nodes.push({p.children_idx[1], node.depth + 1});
      }
    }
    for (auto &node : m_nodes) {
      if (node.is_leaf) {
        node.d_tree.refine();
        // node.d_tree.debugPrint();
      }
    }
    std::cout << "STree refine end\n";
    std::cout << "#BVH nodes: " << m_nodes.size() << std::endl;
  }
  void sample(const Point3f &vertex_pos, const Vector2f &sample,
              SDTreeSampleResult *sds) const {
    auto node = m_nodes[find(vertex_pos)];
    node.d_tree.sample(sample, *sds);
  }
  void update(const Point3f &vertex_pos, const Vector3f &wi, float radiance,
              float pdf, float weight) {
    int node_idx = find(vertex_pos);
    auto &d_tree = m_nodes[node_idx].d_tree;
    d_tree.update(wi, radiance / pdf * weight);
    d_tree.setNumSamples(d_tree.getNumSamples() + weight);
  }

  void debugPrint() const {
    for (int i = 0; i < m_nodes.size(); i++) {
      auto &node = m_nodes[i];
      std::cout << "node\t" << i << ":\n";
      std::cout << "\textend\t" << node.bvh_extend << "\n";
      std::cout << "\taxis_idx\t" << node.axis_idx << "\n";
      std::cout << "\tchildren\t" << node.children_idx[0] << ", "
                << node.children_idx[1] << "\n";
      if (node.is_leaf)
        std::cout << "\t#samples\t" << node.d_tree.getNumSamples() << "\n";
    }
  }

 private:
  void subdivide(int node_index) {
    m_nodes.resize(m_nodes.size() + 2);
    auto &p = m_nodes[node_index];
    auto n_samples = p.d_tree.getNumSamples();

    for (int i = 0; i < 2; i++) {
      int idx = m_nodes.size() - 2 + i;
      p.children_idx[i] = idx;

      m_nodes[idx].axis_idx = (p.axis_idx + 1) % 3;
      m_nodes[idx].d_tree = p.d_tree;
      m_nodes[idx].is_leaf = true;
      m_nodes[idx].d_tree.setNumSamples(n_samples / 2);

      m_nodes[idx].bvh_extend = p.bvh_extend * 0.5f;
    }
    p.is_leaf = false;
    p.d_tree = DTree{};
    /**
     * n_samples = nSamples(leaf)
     * # This modifies the structure. Note the axis.
     * children = splitBinaryTree(leaf)
     * setNSamples(children, n_samples/2)
     */
  }
  int find(const Point3f &vertex_pos) const {
    AABB aabb = m_aabb;  // Shrink during search.
    int cur_node_idx = 0;
    while (!m_nodes[cur_node_idx].is_leaf) {
      auto &cur_node = m_nodes[cur_node_idx];
      Point3f p_min = aabb.pMin;
      Point3f p_max = aabb.pMax;
      Vector3f new_extent = p_max - p_min;
      new_extent[cur_node.axis_idx] *= 0.5f;
      AABB aabb_min{p_min, p_min + new_extent};
      AABB aabb_max{p_max - new_extent, p_max};
      if (aabb_min.Contains(vertex_pos)) {
        aabb = aabb_min;
        cur_node_idx = cur_node.children_idx[0];
      } else {
        aabb = aabb_max;
        cur_node_idx = cur_node.children_idx[1];
      }
    }
    return cur_node_idx;
  }
  std::vector<STreeNode> m_nodes;
  // AABB of root;
  AABB m_aabb;
  int m_sample_balance = 4;
};

void GuidedPathIntegrator::adaptScene(const Scene &scene) {
  assert(!m_sampling_STree && !m_building_STree);
  m_sampling_STree = new STree{scene.boundingBox()};
  m_building_STree = new STree{scene.boundingBox()};
}

void GuidedPathIntegrator::preSampling() {
  // Switch trees.
  // std::swap(m_sampling_STree, m_building_STree);
  STree *t = m_sampling_STree;
  m_sampling_STree = m_building_STree;
  m_building_STree = t;

  m_building_STree->clear();
  m_log << "begin itr " << m_n_iteration << std::endl;
}
void GuidedPathIntegrator::postSampling(int last_itr_spp) {
  // Update is in li, only refine here.
  m_n_iteration++;
  std::cout << "refine building STree " << m_building_STree << std::endl;
  m_building_STree->refine(last_itr_spp);
}

GuidedPathIntegrator::GuidedPathIntegrator(const Json &json)
    : Integrator(json) {
  if (json.contains("max_depth") && json["max_depth"].is_number())
    m_max_depth = json["max_depth"];

  if (json.contains("sample_light") && json["sample_light"].is_boolean())
    m_sample_light = json["sample_light"];

  // Sampling bsdf is replaced by SDTree.

  m_log.open("log.txt", std::ios::out);
}

GuidedPathIntegrator::~GuidedPathIntegrator() {
  m_log.close();
  if (m_sampling_STree) {
    delete m_sampling_STree;
    m_sampling_STree = nullptr;
  }
  if (m_building_STree) {
    delete m_building_STree;
    m_building_STree = nullptr;
  }
}

Spectrum GuidedPathIntegrator::li(Ray &ray, const Scene &scene,
                                  std::shared_ptr<Sampler> sampler) const {
  std::vector<PathVertex> path;
  SDTreeSampleResult sds;
  Vector3f wi;

  Spectrum L{0.f}, beta{1.f};
  int depth = 0;
  // True if current ray is from specular BSDF.
  bool from_specular = true;
  // pdf of sampling this vertex. 1.0 for init ray.
  float vertex_pdf = 1.f;
  do {
    path.push_back({ray.origin, vertex_pdf});
    auto si = scene.rayIntersect(ray);
    if (!si.has_value()) {
      if (!m_sample_light || from_specular)
        // Lights are not sampled until hit,
        // so infinite lights should be added, or they'll never get sampled.
        // Lights are also not sampled if last BSDF is specular,
        // under that condition the direction is fixed.
        for (auto &light : scene.infiniteLights) {
          // Path ends here.
          Spectrum radiance = beta * light->evaluateEmission(ray);
          L += radiance;
          path.back().records.push_back({ray.direction, radiance.Y(), 1.f});
        }
      break;
    }

    // cur_vertex = cur_vertex->next.back();
    auto its = si.value();
    Vector3f wo = -ray.direction;

    if (!m_sample_light || from_specular)  // Direct hit on lights.
      // from_specular also enables direct hit on light.
      if (auto light = its.shape->light; light) {
        Spectrum radiance = beta * light->evaluateEmission(its, wo);
        path.back().records.push_back(
            {ray.direction, radiance.Y(), vertex_pdf});
        L += radiance;
      }
    if (depth++ == m_max_depth) {
      break;
    }
    // Get BSDF.
    auto bsdf = its.shape->material->computeBSDF(its);
    if (!bsdf) {
      // Skip to next intersection in medium rendering.
      ray = Ray(its.position + EPSILON * ray.direction, ray.direction);
      from_specular = true;
      continue;
    }
    // Sample lights as p_{end}.
    // Light sampling starts from this vertex.
    if (m_sample_light) {
      // Infinite lights.
      for (auto light : scene.infiniteLights) {
        auto light_sample = light->sample(its, sampler->next2D());
        Ray shadow_ray(its.position, light_sample.direction, 1e-4f, FLT_MAX);
        // Light connected.
        if (auto occlude = scene.rayIntersect(shadow_ray);
            !occlude.has_value()) {
          Spectrum f = bsdf->f(wo, shadow_ray.direction) *
                       abs(dot(shadow_ray.direction, its.normal));
          float pdf = convertPDF(light_sample, its);

          Spectrum radiance = beta * light_sample.energy * f / pdf;
          path.back().records.push_back(
              {shadow_ray.direction, radiance.Y(), pdf});
          L += radiance;
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

          Spectrum radiance = beta * light_sample.energy * f / pdf;
          path.back().records.push_back(
              {shadow_ray.direction, radiance.Y(), pdf});
          L += radiance;
        }
      }
    }
    // Construct path, accumulate throughput.

    // if (1) {
    if (m_n_iteration == 0) {
      auto bs = bsdf->sample(wo, sampler->next2D());
      wi = bs.wi;
      beta *= bs.weight * abs(dot(bs.wi, its.normal)) / bs.pdf;
      from_specular = bs.type == BSDFType::Specular;
    } else {
      m_sampling_STree->sample(its.position, sampler->next2D(), &sds);
      wi = sds.wi;
      vertex_pdf = sds.pdf;
      beta *= bsdf->f(wo, wi) * abs(dot(wi, its.normal)) / sds.pdf;
      // std::cout << "\n----------------\n";
      // bsdf->f(wo, wi).debugPrint();
      // std::cout << "cos_t " << abs(dot(wi, its.normal)) << std::endl;
      // std::cout << "pdf " << sds.pdf << std::endl;
      // beta.debugPrint();
      // return L;
      from_specular = false;
    }
    ray = Ray{its.position + EPSILON * wi, wi};
    // } while (!beta.isZero());  // Early stop.
  } while (1);

  // Update on-the-fly.
  if (path.size()) {
    /**
     * We got a tree with actual intersection vertex as trunk
     * and light sample directions as branches.
     * We need to count how many times each vertex contributes to the
     * radiance.
     */
    // Accumulate the radiance. Pdf should be multiplied.
    float total_radiance = 0;
    for (int i = path.size() - 1; i >= 0; i--) {
      for (const auto &record : path[i].records) {
        total_radiance += record.radiance;
        m_building_STree->update(path[i].pos, record.dir, record.radiance, 1.f,
                                 1.f);
      }
      if (i) {
        Vector3f dir = path[i].pos - path[i - 1].pos;
        dir = normalize(dir);
        m_building_STree->update(path[i - 1].pos, dir, total_radiance, 1.f,
                                 1.f);
      }
    }
  }
  // return Spectrum{wi[0], wi[1], wi[2]};

  return L;
}

REGISTER_CLASS(GuidedPathIntegrator, "guidedPath")