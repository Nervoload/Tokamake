#include "tokamak/particle_system.hpp"

#include <cmath>

namespace tokamak {

ParticleSystem::ParticleSystem(std::size_t maxParticles) : maxParticles_(maxParticles) {
    positions_.reserve(maxParticles_);
    velocities_.reserve(maxParticles_);
    masses_kg_.reserve(maxParticles_);
    charges_C_.reserve(maxParticles_);
    weights_.reserve(maxParticles_);
    chargeToMass_.reserve(maxParticles_);
    species_.reserve(maxParticles_);
}

bool ParticleSystem::CanInsert(std::size_t count) const {
    if (positions_.size() > maxParticles_) {
        return false;
    }
    return (maxParticles_ - positions_.size()) >= count;
}

bool ParticleSystem::AddParticle(
    const Vec3& position,
    const Vec3& velocity,
    float mass_kg,
    float charge_C,
    ParticleType type,
    float weight) {
    if (!CanInsert(1)) {
        return false;
    }
    if (!position.IsFinite() || !velocity.IsFinite()) {
        return false;
    }
    // Ion-only model in this phase: reject invalid physical particle definitions.
    if (!std::isfinite(mass_kg) || mass_kg <= 0.0f) {
        return false;
    }
    if (!std::isfinite(charge_C) || charge_C == 0.0f) {
        return false;
    }
    if (!std::isfinite(weight) || weight <= 0.0f) {
        return false;
    }

    positions_.push_back(position);
    velocities_.push_back(velocity);
    masses_kg_.push_back(mass_kg);
    charges_C_.push_back(charge_C);
    weights_.push_back(weight);
    chargeToMass_.push_back(charge_C / mass_kg);
    species_.push_back(type);
    return true;
}

bool ParticleSystem::AddParticle(const Vec3& position, const Vec3& velocity, float mass_kg, float charge_C, ParticleType type) {
    return AddParticle(position, velocity, mass_kg, charge_C, type, macroWeight);
}

void ParticleSystem::MarkDead(std::size_t index) {
    species_[index] = ParticleType::Dead;
}

void ParticleSystem::Compact() {
    for (std::size_t i = 0; i < positions_.size();) {
        if (species_[i] == ParticleType::Dead) {
            positions_[i] = positions_.back();
            positions_.pop_back();
            velocities_[i] = velocities_.back();
            velocities_.pop_back();
            masses_kg_[i] = masses_kg_.back();
            masses_kg_.pop_back();
            charges_C_[i] = charges_C_.back();
            charges_C_.pop_back();
            weights_[i] = weights_.back();
            weights_.pop_back();
            chargeToMass_[i] = chargeToMass_.back();
            chargeToMass_.pop_back();
            species_[i] = species_.back();
            species_.pop_back();
        } else {
            ++i;
        }
    }
}

bool ParticleSystem::IsArrayLengthConsistent() const {
    const std::size_t n = positions_.size();
    return velocities_.size() == n && masses_kg_.size() == n && charges_C_.size() == n &&
           weights_.size() == n &&
           chargeToMass_.size() == n && species_.size() == n;
}

bool ParticleSystem::IsFiniteState() const {
    if (!IsArrayLengthConsistent()) {
        return false;
    }

    for (std::size_t i = 0; i < positions_.size(); ++i) {
        if (!positions_[i].IsFinite() || !velocities_[i].IsFinite()) {
            return false;
        }
        if (!std::isfinite(masses_kg_[i]) || !std::isfinite(charges_C_[i]) || !std::isfinite(weights_[i]) ||
            !std::isfinite(chargeToMass_[i]) || weights_[i] <= 0.0f) {
            return false;
        }
    }

    return true;
}

}  // namespace tokamak
