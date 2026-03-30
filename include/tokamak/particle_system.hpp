#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "tokamak/config.hpp"

namespace tokamak {

class ParticleSystem {
public:
    explicit ParticleSystem(std::size_t maxParticles = constants::kDefaultMaxParticles);

    bool CanInsert(std::size_t count) const;
    bool AddParticle(
        const Vec3& position,
        const Vec3& velocity,
        double mass_kg,
        double charge_C,
        ParticleType type,
        double weight);
    bool AddParticle(const Vec3& position, const Vec3& velocity, double mass_kg, double charge_C, ParticleType type);
    void MarkDead(std::size_t index);
    void RestoreSpecies(std::size_t index, ParticleType type);
    void Compact();

    std::size_t Size() const { return positions_.size(); }
    std::size_t MaxParticles() const { return maxParticles_; }

    bool IsArrayLengthConsistent() const;
    bool IsFiniteState() const;

    const std::vector<Vec3>& Positions() const { return positions_; }
    const std::vector<Vec3>& Velocities() const { return velocities_; }
    const std::vector<double>& Masses() const { return masses_kg_; }
    const std::vector<double>& Charges() const { return charges_C_; }
    const std::vector<double>& Weights() const { return weights_; }
    const std::vector<double>& ChargeToMass() const { return chargeToMass_; }
    const std::vector<uint64_t>& Ids() const { return ids_; }
    const std::vector<ParticleType>& Species() const { return species_; }

    std::vector<Vec3>& MutablePositions() { return positions_; }
    std::vector<Vec3>& MutableVelocities() { return velocities_; }
    std::vector<double>& MutableMasses() { return masses_kg_; }
    std::vector<double>& MutableCharges() { return charges_C_; }
    std::vector<double>& MutableWeights() { return weights_; }
    std::vector<double>& MutableChargeToMass() { return chargeToMass_; }
    std::vector<uint64_t>& MutableIds() { return ids_; }
    std::vector<ParticleType>& MutableSpecies() { return species_; }

    ParticleType SpeciesAt(std::size_t index) const { return species_.at(index); }

    double macroWeight = 1.0e12;
    uint64_t fusionCountTotal = 0;

private:
    std::size_t maxParticles_;
    std::size_t deadParticleCount_ = 0;
    std::vector<Vec3> positions_;
    std::vector<Vec3> velocities_;
    std::vector<double> masses_kg_;
    std::vector<double> charges_C_;
    std::vector<double> weights_;
    std::vector<double> chargeToMass_;
    std::vector<uint64_t> ids_;
    std::vector<ParticleType> species_;
    uint64_t nextParticleId_ = 1;
};

}  // namespace tokamak
