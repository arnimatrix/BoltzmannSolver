#pragma once

#include <marty.h>

namespace mty {

    class QuantumNumberData;

    class ProcessData {

    private:

        struct Process {
            std::string name;
            std::vector<mty::Insertion> inParticles;
            std::vector<mty::Insertion> outParticles;
            std::map<std::string, int> qNumbers;
        };

    public:

        ProcessData(mty::Model const &t_model);

        void addProcess(
                std::string                 const &name,
                std::vector<mty::Insertion> const &insertions,
                QuantumNumberData           const &qData
                );

        std::unique_ptr<JSON::Node> generateJSONParticleList(
                std::string                 const &listName,
                std::vector<mty::Insertion> const &particles
                ) const;

        std::unique_ptr<JSON::Node> generateJSONQNumberList(
                std::string                const &listName,
                std::map<std::string, int> const &qNumbers
                ) const;

        std::unique_ptr<JSON::Node> generateJSONData() const;

    private:

        void processQuantumNumbers(
                Process                 &process,
                QuantumNumberData const &qData
                ) const;

    private:

        mty::Model const *m_model;
        std::vector<Process> processes;
    };

    void saveParticleData(
            std::string       const &fileName,
            QuantumNumberData const &qData,
            ProcessData       const &pData
            );
}
