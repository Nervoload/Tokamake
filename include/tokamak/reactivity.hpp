#pragma once

#include <vector>

#include "tokamak/config.hpp"

namespace tokamak {

struct SigmaETablePoint {
    float energy_keV = 0.0f;
    float sigma_m2 = 0.0f;
};

using SigmaETable = std::vector<SigmaETablePoint>;

SigmaETable MakeDefaultDtSigmaETable();
bool IsValidSigmaETable(const SigmaETable& table);

double EvaluateSigmaFromTable(double energy_keV, const SigmaETable& table);
double EvaluateDtSigma_m2(double energy_keV, FusionReactivityModelKind kind);

}  // namespace tokamak
