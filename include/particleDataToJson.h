#pragma once

#include <string>
#include <vector>
#include <memory>

namespace JSON {
  class Object;
  class Node;
  class List;
} // namespace JSON

namespace mty {

  class Particle;
  class Model;
  class QuantumNumber;

  std::unique_ptr<JSON::Object> quantumNumberToJson(
      mty::QuantumNumber const &qNumber,
      int                       value
      );

  std::unique_ptr<JSON::Node> particleDataToJson(
      std::string                     const &specifier,
      mty::Particle                   const &particle,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      );

  std::unique_ptr<JSON::Node> particleToJson(
      std::string                     const &specifier,
      mty::Particle                   const &particle,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      );

  std::unique_ptr<JSON::List> particleListToJson(
      std::vector<mty::Particle>      const &particles,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      );

  std::unique_ptr<JSON::Node> modelToJsonData(
      mty::Model const &model
      );

  void saveModelDatatoJson(
      std::string const &fileName,
      mty::Model  const &model
      );

} // namespace mty
