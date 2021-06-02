#include "particleData.h"
#include <fstream>
#include <algorithm>
#include "marty/jsonLoader.h"
#include "marty/jsonObject.h"

namespace mty {

    ParticleData ParticleData::readFromJson(std::string const &fileName)
    {
        using namespace JSON;

        ParticleData res;
        std::unique_ptr<Node> data = Reader::loadFromFile(fileName);
        auto particleList = Parser::parseNode(data.get(), "particles", true);
        for (const auto &particle : *particleList) {
            std::string name = Parser::parseArgument<std::string>(
                    particle, "name", true).value();
            res.addParticle(name);
            auto qNumbers = Parser::parseNode(particle, "qnumbers", true);
            for (const auto &qN : *qNumbers) {
                std::string qName = qN->getSpecifier();
                int         value = Parser::interpretObject<int>(qN.get());
                res.addQuantumNumber(qName);
                res(qName, name) = value;
            }
        }

        return res;
    }

    void ParticleData::addElement(
            std::string        const &element,
            std::vector<std::string> &vec
            )
    {
        auto pos = std::find(vec.begin(), vec.end(), element);
        if (pos == vec.end()) {
            vec.push_back(element);
        }
    }

    std::ostream &operator<<(
            std::ostream       &out,
            ParticleData const &data
            )
    {
        out << "ParticleData object\n";
        out << "Quantum numbers:\n";
        for (const auto &qNumber : data.getQuantumNumbers()) {
            out << "\t-> " << qNumber << '\n';
        }
        out << "Particles:\n";
        for (const auto &particle : data.getParticleNames()) {
            out << "\t-> " << particle << '\n';
        }
        out << "Map:\n";
        for (const auto &[qNumberData, value] : data.data) {
            out << "\t-> " << qNumberData.qNumber << ", "
                << qNumberData.particle << " = " << value << '\n';
        }
        out.flush();
        return out;
    }

} // namespace mty
