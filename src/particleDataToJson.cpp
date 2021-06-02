#include "particleDataToJson.h"
#include <marty.h>
#include <marty/jsonLoader.h>
#include <marty/jsonObject.h>

namespace mty {

  std::unique_ptr<JSON::Object> quantumNumberToJson(
      mty::QuantumNumber const &qNumber,
      int                       value
      )
  {
      return JSON::Leaf<int>::make(qNumber.getName(), value);
  }

  std::unique_ptr<JSON::Node> particleDataToJson(
      std::string                     const &specifier,
      mty::Particle                   const &particle,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      )
  {
      std::unique_ptr<JSON::Node> qList = JSON::Node::make(specifier);
      for (const auto &qNumber : quantumNumbers) {
          std::unique_ptr<JSON::Object> node = quantumNumberToJson(
              qNumber, particle->getQuantumNumber(&qNumber)
              );
          qList->addChild(node);
      }
      return qList;
  }
  
  std::unique_ptr<JSON::Node> particleToJson(
      std::string                     const &specifier,
      mty::Particle                   const &particle,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      )
  {
      std::unique_ptr<JSON::Node> node = JSON::Node::make(specifier);
      auto nameLeaf = JSON::Leaf<std::string>::make("name", particle->getName());
      auto qnumberList = particleDataToJson(
          "qnumbers", particle, quantumNumbers);
      node->addChild(nameLeaf);
      node->addChild(qnumberList);
      return node;
  }

  std::unique_ptr<JSON::List> particleListToJson(
      std::vector<mty::Particle>      const &particles,
      std::vector<mty::QuantumNumber> const &quantumNumbers
      )
  {
      constexpr auto specifier = "particles";
      std::unique_ptr<JSON::List> li = JSON::List::make(specifier);
      for (const auto &part : particles) {
          auto particleNode = particleToJson(specifier, part, quantumNumbers);
          li->addChild(particleNode);
      }
      return li;
  }

  std::unique_ptr<JSON::Node> modelToJsonData(
      mty::Model const &model
      )
  {
      std::unique_ptr<JSON::Node> root = JSON::Node::make("root");
      auto particleList = particleListToJson(
          model.getParticles(),
          model.getQuantumNumbers()
          );
      root->addChild(particleList);
      return root;
  }

  void saveModelDatatoJson(
      std::string const &fileName,
      mty::Model  const &model
      )
  {
      JSON::Reader::saveToFile(fileName, modelToJsonData(model));
  }

} // namespace mty
