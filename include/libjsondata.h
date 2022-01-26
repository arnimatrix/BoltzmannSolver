#pragma once

#include <functional>
#include <string>
#include <vector>
#include <map>
#include "marty/jsonLoader.h"

namespace mty::lib {

    struct Particle {
        std::string name;
        bool        particle { true };

        static Particle buildFromSpecifier(std::string spec) {
            size_t const sz = spec.size();
            if (sz > 2 && spec[sz-2] == '#' && spec[sz-1] == 'c') {
                spec.erase(sz-2, 2);
                return {spec, false};
            }
            return {spec, true};
        }

        bool operator==(Particle const &other) const {
            return name == other.name && particle == other.particle;
        }
    };

    struct Process {
        std::string name;
        std::vector<Particle> inParticles;
        std::vector<Particle> outParticles;
        std::map<std::string, int> qNumbers;
    };

    struct Rational {
        int num;
        int denom;

        double getValue() const { return num * 1. / denom; }
    };

    struct QuantumNumber {
        std::string name;
        int         factor;
        std::map<std::string, int> particles;

        Rational getQuantumNumberValue(
            std::string const &particleName,
            bool               isParticle = true
            );

        Rational getQuantumNumberValue(
            Particle const &particle
            );
    };

    std::ostream &operator<<(std::ostream &out, Process const &process);
    std::ostream &operator<<(std::ostream &out, QuantumNumber const &qNumber);

    class ParticleData {

        public:

        std::vector<std::string> const &getParticleNames() const {
            return particleNames;
        }

        std::vector<QuantumNumber> const &getQuantumNumbers() const {
            return qNumbers;
        }

        Rational getQuantumNumberValue(
            std::string const &quantumNumber,
            std::string const &particleName,
            bool               isParticle = true
            );

        Rational getQuantumNumberValue(
            std::string const &quantumNumber,
            Particle    const &particle
            );

        std::vector<Process> getProcesses(
                std::function<bool(Process const&)> selection
                ) const;

        std::vector<Process> getProcessesFromIncoming(
                std::string const &p1,
                std::string const &p2 = ""
                ) const;

        std::vector<Process> getProcessesFromQNumberConservation(
                std::string const &qNumber
                ) const;

        std::vector<Process> getProcessesFromQNumberViolation(
                std::string const &qNumber
                ) const;

        void loadFile(std::string const &nameFile);

        private:

        void clear();

        void load(std::unique_ptr<JSON::Node> const &tree);

        void loadQNumberNode(JSON::Node const *node);

        void loadProcessNode(JSON::Node const *node);

        bool isConjugated(std::string &name) const;

        std::vector<Particle> loadParticles(JSON::Node const *node) const;

        void addParticles(std::vector<Particle> const &particles);

        void loadProcessQNumber(
                JSON::Node const *node,
                Process          &process
                ) const;

        friend
            std::ostream &operator<<(
                    std::ostream       &out,
                    ParticleData const &data
                    );

        private:

        std::vector<std::string> particleNames;

        std::vector<QuantumNumber> qNumbers;

        std::vector<Process> processes;
    };

}
