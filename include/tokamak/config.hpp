#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

#include "tokamak/types.hpp"

namespace tokamak {

namespace constants {

constexpr float kPi = 3.14159265359f;
constexpr float kMu0 = 1.25663706e-6f;
constexpr float kElementaryCharge_C = 1.60217663e-19f;
constexpr float kMassDeuterium_kg = 3.3435e-27f;
constexpr float kMassTritium_kg = 5.0082e-27f;
constexpr float kMassHelium4_kg = 6.6464e-27f;
constexpr float kBoltzmann_JPerK = 1.380649e-23f;
constexpr std::size_t kDefaultMaxParticles = 2500000;

}  // namespace constants

enum class Scenario : uint8_t {
    ColdVacuum = 0,
    NbiIgnition = 1,
    MagneticFailure = 2,
};

enum class ParticleType : uint8_t {
    Deuterium = 0,
    Tritium = 1,
    Helium = 2,
    Dead = 3,
};

enum class PlasmaCurrentProfileKind : uint8_t {
    Uniform = 0,
    Parabolic = 1,
    CustomTable = 2,
};

enum class ElectricFieldMode : uint8_t {
    Placeholder = 0,
    Electrostatic = 1,
};

enum class ElectrostaticBoundaryCondition : uint8_t {
    DirichletZero = 0,
    NeumannZeroGradient = 1,
};

enum class ChargeAssignmentScheme : uint8_t {
    NGP = 0,
    CIC = 1,
};

enum class FusionReactivityModelKind : uint8_t {
    SigmaE_Table = 0,
};

enum class WallBoundaryMode : uint8_t {
    Reflect = 0,
    Absorb = 1,
    Recycle = 2,
};

struct CurrentProfilePoint {
    float normalizedMinorRadius = 0.0f;
    float enclosedCurrentFraction = 0.0f;
};

inline const char* ScenarioName(Scenario scenario) {
    switch (scenario) {
        case Scenario::ColdVacuum:
            return "COLD_VACUUM";
        case Scenario::NbiIgnition:
            return "NBI_IGNITION";
        case Scenario::MagneticFailure:
            return "MAGNETIC_FAILURE";
    }
    return "UNKNOWN";
}

inline bool ParseScenario(std::string_view text, Scenario* outScenario) {
    if (text == "cold" || text == "COLD_VACUUM") {
        *outScenario = Scenario::ColdVacuum;
        return true;
    }
    if (text == "ignition" || text == "NBI_IGNITION") {
        *outScenario = Scenario::NbiIgnition;
        return true;
    }
    if (text == "failure" || text == "MAGNETIC_FAILURE") {
        *outScenario = Scenario::MagneticFailure;
        return true;
    }
    return false;
}

inline const char* PlasmaCurrentProfileKindName(PlasmaCurrentProfileKind kind) {
    switch (kind) {
        case PlasmaCurrentProfileKind::Uniform:
            return "uniform";
        case PlasmaCurrentProfileKind::Parabolic:
            return "parabolic";
        case PlasmaCurrentProfileKind::CustomTable:
            return "custom";
    }
    return "unknown";
}

inline bool ParsePlasmaCurrentProfileKind(std::string_view text, PlasmaCurrentProfileKind* outKind) {
    if (text == "uniform" || text == "UNIFORM") {
        *outKind = PlasmaCurrentProfileKind::Uniform;
        return true;
    }
    if (text == "parabolic" || text == "PARABOLIC") {
        *outKind = PlasmaCurrentProfileKind::Parabolic;
        return true;
    }
    if (text == "custom" || text == "CUSTOM") {
        *outKind = PlasmaCurrentProfileKind::CustomTable;
        return true;
    }
    return false;
}

inline const char* ElectricFieldModeName(ElectricFieldMode mode) {
    switch (mode) {
        case ElectricFieldMode::Placeholder:
            return "placeholder";
        case ElectricFieldMode::Electrostatic:
            return "electrostatic";
    }
    return "unknown";
}

inline bool ParseElectricFieldMode(std::string_view text, ElectricFieldMode* outMode) {
    if (text == "placeholder" || text == "PLACEHOLDER") {
        *outMode = ElectricFieldMode::Placeholder;
        return true;
    }
    if (text == "electrostatic" || text == "ELECTROSTATIC") {
        *outMode = ElectricFieldMode::Electrostatic;
        return true;
    }
    return false;
}

inline const char* ElectrostaticBoundaryConditionName(ElectrostaticBoundaryCondition boundaryCondition) {
    switch (boundaryCondition) {
        case ElectrostaticBoundaryCondition::DirichletZero:
            return "dirichlet0";
        case ElectrostaticBoundaryCondition::NeumannZeroGradient:
            return "neumann0";
    }
    return "unknown";
}

inline bool ParseElectrostaticBoundaryCondition(std::string_view text, ElectrostaticBoundaryCondition* outBoundaryCondition) {
    if (text == "dirichlet0" || text == "DIRICHLET0") {
        *outBoundaryCondition = ElectrostaticBoundaryCondition::DirichletZero;
        return true;
    }
    if (text == "neumann0" || text == "NEUMANN0") {
        *outBoundaryCondition = ElectrostaticBoundaryCondition::NeumannZeroGradient;
        return true;
    }
    return false;
}

inline const char* ChargeAssignmentSchemeName(ChargeAssignmentScheme scheme) {
    switch (scheme) {
        case ChargeAssignmentScheme::NGP:
            return "ngp";
        case ChargeAssignmentScheme::CIC:
            return "cic";
    }
    return "unknown";
}

inline bool ParseChargeAssignmentScheme(std::string_view text, ChargeAssignmentScheme* outScheme) {
    if (text == "ngp" || text == "NGP") {
        *outScheme = ChargeAssignmentScheme::NGP;
        return true;
    }
    if (text == "cic" || text == "CIC") {
        *outScheme = ChargeAssignmentScheme::CIC;
        return true;
    }
    return false;
}

inline const char* FusionReactivityModelKindName(FusionReactivityModelKind kind) {
    switch (kind) {
        case FusionReactivityModelKind::SigmaE_Table:
            return "sigmae-table";
    }
    return "unknown";
}

inline bool ParseFusionReactivityModelKind(std::string_view text, FusionReactivityModelKind* outKind) {
    if (text == "sigmae-table" || text == "SIGMAE_TABLE" || text == "SIGMAE-TABLE") {
        *outKind = FusionReactivityModelKind::SigmaE_Table;
        return true;
    }
    return false;
}

inline const char* WallBoundaryModeName(WallBoundaryMode mode) {
    switch (mode) {
        case WallBoundaryMode::Reflect:
            return "reflect";
        case WallBoundaryMode::Absorb:
            return "absorb";
        case WallBoundaryMode::Recycle:
            return "recycle";
    }
    return "unknown";
}

inline bool ParseWallBoundaryMode(std::string_view text, WallBoundaryMode* outMode) {
    if (text == "reflect" || text == "REFLECT") {
        *outMode = WallBoundaryMode::Reflect;
        return true;
    }
    if (text == "absorb" || text == "ABSORB") {
        *outMode = WallBoundaryMode::Absorb;
        return true;
    }
    if (text == "recycle" || text == "RECYCLE") {
        *outMode = WallBoundaryMode::Recycle;
        return true;
    }
    return false;
}

struct TokamakConfig {
    float majorRadius_m = 2.0f;
    float minorRadius_m = 0.5f;
    float toroidalCurrent_A = 15.0e6f;
    int toroidalCoilTurns = 18;
    float plasmaCurrent_A = 2.0e6f;
};

struct NBIConfig {
    bool isActive = true;
    float beamEnergy_keV = 100.0f;
    int particlesPerStep = 50;
    Vec3 injectorPos = Vec3(2.5f, 0.0f, 0.0f);
    Vec3 injectionNormal = Vec3(-1.0f, 0.2f, 0.0f).Normalized();
};

struct RunConfig {
    Scenario scenario = Scenario::NbiIgnition;
    float timeStep_s = 1.0e-7f;
    int totalSteps = 10000;
    int telemetryEveryNSteps = 100;
    std::optional<uint32_t> seed;
    std::size_t particleCap = constants::kDefaultMaxParticles;

    PlasmaCurrentProfileKind plasmaCurrentProfileKind = PlasmaCurrentProfileKind::Uniform;
    float plasmaCurrentAxisEpsilon_m = 1.0e-4f;
    float plasmaCurrentCustomAxisBlendRadius_m = 1.0e-3f;
    std::vector<CurrentProfilePoint> plasmaCurrentCustomTable;

    std::size_t magneticFieldRadialBinCount = 32;
    double magneticFieldDtSafetyFraction = 0.05;

    ElectricFieldMode electricFieldMode = ElectricFieldMode::Placeholder;
    ElectrostaticBoundaryCondition electrostaticBoundaryCondition = ElectrostaticBoundaryCondition::DirichletZero;
    ChargeAssignmentScheme chargeAssignmentScheme = ChargeAssignmentScheme::CIC;
    std::size_t electrostaticGridBinCount = 32;
    double electrostaticSolverTolerance = 1.0e-6;
    uint32_t electrostaticSolverMaxIterations = 5000;
    double electrostaticSorOmega = 1.7;
    double electrostaticNeutralizingBackgroundFraction = 1.0;

    FusionReactivityModelKind fusionReactivityModelKind = FusionReactivityModelKind::SigmaE_Table;
    double fusionCrossSectionScale = 1.0;
    double fusionProbabilityClamp = 0.95;
    double fusionMinEnergy_keV = 0.0;
    WallBoundaryMode wallBoundaryMode = WallBoundaryMode::Reflect;
    double recycleFraction = 0.0;
    std::size_t fusionDiagnosticsRadialBins = 32;
};

}  // namespace tokamak
