#include "tokamak/reactivity.hpp"

#include <algorithm>
#include <cmath>

namespace tokamak {
namespace {

constexpr double kBarnToSquareMeter = 1.0e-28;

const SigmaETable& DefaultDtSigmaETableRef() {
    static const SigmaETable kTable = MakeDefaultDtSigmaETable();
    return kTable;
}

}  // namespace

SigmaETable MakeDefaultDtSigmaETable() {
    return SigmaETable{
        SigmaETablePoint{0.0f, 0.0f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{5.0f, 0.01f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{10.0f, 0.05f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{20.0f, 0.20f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{30.0f, 0.45f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{50.0f, 1.20f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{80.0f, 2.20f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{120.0f, 3.00f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{200.0f, 3.60f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{300.0f, 4.00f * static_cast<float>(kBarnToSquareMeter)},
        SigmaETablePoint{500.0f, 4.30f * static_cast<float>(kBarnToSquareMeter)},
    };
}

bool IsValidSigmaETable(const SigmaETable& table) {
    if (table.size() < 2) {
        return false;
    }
    for (std::size_t i = 0; i < table.size(); ++i) {
        const auto& point = table[i];
        if (!std::isfinite(point.energy_keV) || !std::isfinite(point.sigma_m2) || point.energy_keV < 0.0f ||
            point.sigma_m2 < 0.0f) {
            return false;
        }
        if (i > 0 && table[i - 1].energy_keV >= point.energy_keV) {
            return false;
        }
    }
    return true;
}

double EvaluateSigmaFromTable(double energy_keV, const SigmaETable& table) {
    if (!IsValidSigmaETable(table) || !std::isfinite(energy_keV)) {
        return 0.0;
    }

    if (energy_keV <= static_cast<double>(table.front().energy_keV)) {
        return static_cast<double>(table.front().sigma_m2);
    }
    if (energy_keV >= static_cast<double>(table.back().energy_keV)) {
        return static_cast<double>(table.back().sigma_m2);
    }

    for (std::size_t i = 1; i < table.size(); ++i) {
        const double e0 = static_cast<double>(table[i - 1].energy_keV);
        const double e1 = static_cast<double>(table[i].energy_keV);
        if (energy_keV > e1) {
            continue;
        }

        const double s0 = static_cast<double>(table[i - 1].sigma_m2);
        const double s1 = static_cast<double>(table[i].sigma_m2);
        const double t = (energy_keV - e0) / (e1 - e0);
        return s0 + (t * (s1 - s0));
    }

    return static_cast<double>(table.back().sigma_m2);
}

double EvaluateDtSigma_m2(double energy_keV, FusionReactivityModelKind kind) {
    switch (kind) {
        case FusionReactivityModelKind::SigmaE_Table:
            return EvaluateSigmaFromTable(energy_keV, DefaultDtSigmaETableRef());
    }
    return 0.0;
}

}  // namespace tokamak
