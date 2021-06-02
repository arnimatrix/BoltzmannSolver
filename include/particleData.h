#pragma once

#include <map>
#include <vector>
#include <string>
#include <string_view>

namespace mty {

    struct QuantumNumberData {
        std::string qNumber;
        std::string particle;

        bool operator==(QuantumNumberData const &other) const {
            return qNumber == other.qNumber
                && particle == other.particle;
        }

        bool operator<(QuantumNumberData const &other) const {
            return qNumber < other.qNumber 
                || (qNumber == other.qNumber && particle < other.particle);
        }
    };

    class ParticleData {

    public:

        std::vector<std::string> const &getParticleNames() const {
            return particleNames;
        }

        std::vector<std::string> const &getQuantumNumbers() const {
            return quantumNumbers;
        }

        bool  empty() const { return data.empty(); }
        size_t size() const { return data.size();  }

        auto begin()       { return data.begin(); }
        auto begin() const { return data.begin(); }
        auto end()         { return data.end();   }
        auto end()   const { return data.end();   }

        auto find(
                std::string const &quantumNumber,
                std::string const &particle
                ) 
        {
            return data.find({quantumNumber, particle});
        }

        auto find(
                std::string const &quantumNumber,
                std::string const &particle
                ) const 
        {
            return data.find({quantumNumber, particle});
        }

        int &operator()(
                std::string const &quantumNumber,
                std::string const &particle
                )
        {
            return data[{quantumNumber, particle}];
        }

        template<class ...Args>
        auto insert(Args &&...args) {
            return data.insert(std::forward<Args>(args)...);
        }

        static ParticleData readFromJson(std::string const &fileName);

        friend std::ostream &operator<<(
                std::ostream       &out,
                ParticleData const &data
                );

    private:

        static void addElement(
                std::string        const &el,
                std::vector<std::string> &vec
                );
        void addQuantumNumber(std::string const &qNumber) {
            addElement(qNumber, quantumNumbers);
        }
        void addParticle     (std::string const &particle) {
            addElement(particle, particleNames);
        }

    private:

        std::map<QuantumNumberData, int> data;

        std::vector<std::string> quantumNumbers;
        std::vector<std::string> particleNames;
    };

} // namespace mty
