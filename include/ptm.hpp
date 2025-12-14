#pragma once

#include <limits>
#include "types.hpp"

namespace wasabi::ptm {
enum class PtmTarget { Residue, ProteinNTerm, ProteinCTerm };
enum class PtmType { Fixed, Variable };

// Number of PTMs in the common list. The type can be adjusted in the future releases if needed.
using common_ptms_size_t = uint8_t;
constexpr common_ptms_size_t common_ptms_size_k = 3;
constexpr common_ptms_size_t no_op_ptm_idx = std::numeric_limits<common_ptms_size_t>::max();

struct PTM {
    const char* name;
    size_t name_size;
    PtmTarget targetType;
    char targetResidue;
    Composition delta;
    PtmType type;
    common_ptms_size_t common_index;
};

// Common PTM constants
constexpr PTM CARBAMIDOMETHYL_C = {
    "Carbamidomethyl", 15,PtmTarget::Residue, 'C', Composition{2,3,1,1,0,0, 0}, PtmType::Fixed, 0
};

constexpr PTM OXIDATION_M = {
    "Oxidation", 9, PtmTarget::Residue, 'M', Composition{0,0,0,1,0,0,0}, PtmType::Variable, 1
};

constexpr PTM ACETYL_NTERM = {
    "Acetyl", 6, PtmTarget::ProteinNTerm, '\0', Composition{2,2,0,1,0,0,0}, PtmType::Variable, 2
};

constexpr PTM PHOSPHORYLATION_M = {
    "Phosphorylation", 15, PtmTarget::Residue, 'M', Composition{0,1,0,3,1,0, 0}, PtmType::Variable, no_op_ptm_idx
};

// A ready-to-use global PTM list
constexpr std::array<PTM, common_ptms_size_k> common_ptms_k{
    CARBAMIDOMETHYL_C,
    OXIDATION_M,
    ACETYL_NTERM
};
// Index in this table identifies PTM. Not registered PTMs are virtually non-existent.

// Function should be updated with each addition of new PTMs.
inline std::pair<size_t, common_ptms_size_t> decode_modification(const std::string& mod) {
    common_ptms_size_t index;
    switch (mod[0]) {
        case 'C':
            index = CARBAMIDOMETHYL_C.common_index;
            break;
        case 'O':
            index = OXIDATION_M.common_index;
            break;
        case 'A':
            index = ACETYL_NTERM.common_index;
            break;
        default:
            throw std::runtime_error("Invalid Modification: " + mod);
    }
    if (mod.size() == common_ptms_k[index].name_size) {
        return std::make_pair(1, index);
    }

    size_t times = 0;
    std::string number_rev;
    for (int i = mod.size() - 1; i >= 0; --i) {
        if (mod[i] == 'x') {
            break;
        }
        number_rev.push_back(mod[i]);
    }
    std::ranges::reverse(number_rev);
    times = static_cast<size_t>(std::stoull(number_rev));

    return std::make_pair(times, index);
}

}