#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Camera/Pinhole.h>
#include <FunctionLayer/Integrator/Integrator.h>
#include <FunctionLayer/Sampler/Sampler.h>
#include <FunctionLayer/Scene/Scene.h>
#include <FunctionLayer/Texture/Mipmap.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/FileUtil.h>
#include <ResourceLayer/Image.h>
#include <ResourceLayer/JsonUtil.h>
#include <stdio.h>

#include <chrono>
#include <fstream>
#include <regex>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(float percentage) {
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

void save() {}

void render(std::shared_ptr<Integrator> integrator,
            std::shared_ptr<Sampler> sampler, std::shared_ptr<Camera> camera,
            std::shared_ptr<Scene> scene, int spp) {
  integrator->preSampling();

  // int spp = sampler->xSamples * sampler->ySamples;
  int width = camera->film->size[0], height = camera->film->size[1];

  auto start = std::chrono::system_clock::now();
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      Vector2f NDC{(float)x / width, (float)y / height};
      Spectrum li(.0f);
      for (int i = 0; i < spp; ++i) {
        Ray ray = camera->sampleRayDifferentials(
            CameraSample{sampler->next2D()}, NDC);
        li += integrator->li(ray, *scene, sampler);
      }
      camera->film->deposit({x, y}, li / spp);

      int finished = x + y * width;
      if (finished % 5 == 0) {
        printProgress((float)finished / (height * width));
      }
    }
  }
  printProgress(1.f);
  auto end = std::chrono::system_clock::now();
  printf("\nRendering costs %.2fs\n",
         (std::chrono::duration_cast<std::chrono::milliseconds>(end - start))
                 .count() /
             1000.f);

  printf("tracing record:\n");
  for (auto& kv : integrator->record.entries) {
    printf("%s:\t%d\n", kv.first.c_str(), kv.second);
  }

  std::cout << "post sampling\n";
  integrator->postSampling(spp);
}

std::string getSuffix(const std::string& outputName, int n_itr, int spp) {
  char suffix_buf[20];
  std::sprintf(suffix_buf, "%04d", n_itr);
  std::string suffix = std::string(suffix_buf) + ".spp" + std::to_string(spp);
  size_t dot_pos = outputName.find_last_of(".");
  std::string ret =
      outputName.substr(0, dot_pos + 1) + suffix + outputName.substr(dot_pos);
  return ret;
}

int main(int argc, char** argv) {
  const std::string sceneDir = std::string(argv[1]);
  FileUtil::setWorkingDirectory(sceneDir);
  std::string sceneJsonPath = FileUtil::getFullPath("scene.json");
  std::ifstream fstm(sceneJsonPath);
  Json json = Json::parse(fstm);
  auto camera = Factory::construct_class<Camera>(json["camera"]);
  auto scene = std::make_shared<Scene>(json["scene"]);
  auto integrator = Factory::construct_class<Integrator>(json["integrator"]);
  auto sampler = Factory::construct_class<Sampler>(json["sampler"]);

  integrator->adaptScene(*scene);
  int n_itr = fetchRequired<int>(json["integrator"], "iteration");
  int spp = sampler->xSamples * sampler->ySamples;
  for (int i = 1; i <= n_itr; i++) {
    std::cout << "Iteration " << i << "/" << n_itr << ", spp " << spp
              << std::endl;
    render(integrator, sampler, camera, scene, spp);
    std::cout << "Saving to file." << std::endl;
    //* 目前支持输出为png/hdr两种格式
    std::string outputName =
        fetchRequired<std::string>(json["output"], "filename");
    outputName = getSuffix(outputName, i, spp);

    std::cout << "Saving to " << outputName << std::endl;
    if (std::regex_match(outputName, std::regex("(.*)(\\.png)"))) {
      camera->film->savePNG(outputName.c_str());
    } else if (std::regex_match(outputName, std::regex("(.*)(\\.hdr)"))) {
      camera->film->saveHDR(outputName.c_str());
    } else {
      std::cout << "Only support output as PNG/HDR\n";
    }
    std::cout << "Iteration done\n";
    spp *= 2;
  }

  return 0;
}
