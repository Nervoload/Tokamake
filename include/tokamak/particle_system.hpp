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
    bool AddParticle(const Vec3& position, const Vec3& velocity, float mass_kg, float charge_C, ParticleType type);
    void MarkDead(std::size_t index);
    void Compact();

    std::size_t Size() const { return positions_.size(); }
    std::size_t MaxParticles() const { return maxParticles_; }

    bool IsArrayLengthConsistent() const;
    bool IsFiniteState() const;

    const std::vector<Vec3>& Positions() const { return positions_; }
    const std::vector<Vec3>& Velocities() const { return velocities_; }
    const std::vector<float>& Masses() const { return masses_kg_; }
    const std::vector<float>& Charges() const { return charges_C_; }
    const std::vector<float>& ChargeToMass() const { return chargeToMass_; }
    const std::vector<ParticleType>& Species() const { return species_; }

    std::vector<Vec3>& MutablePositions() { return positions_; }
    std::vector<Vec3>& MutableVelocities() { return velocities_; }
    std::vector<float>& MutableMasses() { return masses_kg_; }
    std::vector<float>& MutableCharges() { return charges_C_; }
    std::vector<float>& MutableChargeToMass() { return chargeToMass_; }
    std::vector<ParticleType>& MutableSpecies() { return species_; }

    ParticleType SpeciesAt(std::size_t index) const { return species_.at(index); }

    float macroWeight = 1.0e12f;
    uint64_t fusionCountTotal = 0;

private:
    std::size_t maxParticles_;
    std::vector<Vec3> positions_;
    std::vector<Vec3> velocities_;
    std::vector<float> masses_kg_;
    std::vector<float> charges_C_;
    std::vector<float> chargeToMass_;
    std::vector<ParticleType> species_;
};

}  // namespace tokamak
