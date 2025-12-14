#pragma once

#include <array>

#include "types.hpp"

namespace wasabi::aa {

constexpr int kAAOffset = 'A';
constexpr std::size_t kAALength = 'Z' - 'A' + 1;

// Amino acid residues without water
// {C, H, N, O, S, P, Se}
constexpr std::array<Composition, kAALength> kAminoAcidCompositions{{
    /* A */ Composition(3,5,1,1,0,0,0),
    /* B */ Composition(4,5,1,3,0,0,0),
    /* C */ Composition(3,5,1,1,1,0,0),
    /* D */ Composition(4,5,1,3,0,0,0),
    /* E */ Composition(5,7,1,3,0,0,0),
    /* F */ Composition(9,9,1,1,0,0,0),
    /* G */ Composition(2,3,1,1,0,0,0),
    /* H */ Composition(6,7,3,1,0,0,0),
    /* I */ Composition(6,11,1,1,0,0,0),
    /* J */ Composition(6,11,1,1,0,0,0),
    /* K */ Composition(6,12,2,1,0,0,0),
    /* L */ Composition(6,11,1,1,0,0,0),
    /* M */ Composition(5,9,1,1,1,0,0),
    /* N */ Composition(4,6,2,2,0,0,0),
    /* O */ Composition(12,19,3,2,0,0,0),
    /* P */ Composition(5,7,1,1,0,0,0),
    /* Q */ Composition(5,8,2,2,0,0,0),
    /* R */ Composition(6,12,4,1,0,0,0),
    /* S */ Composition(3,5,1,2,0,0,0),
    /* T */ Composition(4,7,1,2,0,0,0),
    /* U */ Composition(3,5,1,1,0,0,1),
    /* V */ Composition(5,9,1,1,0,0,0),
    /* W */ Composition(11,10,2,1,0,0,0),
    /* X */ Composition(5,9,1,1,0,0,0),
    /* Y */ Composition(9,9,1,2,0,0,0),
    /* Z */ Composition(5,7,1,3,0,0,0)
}};

constexpr const Composition& GetAminoAcidComposition(const char aa) noexcept {
    return (aa >= 'A' && aa <= 'Z')
        ? kAminoAcidCompositions[aa - kAAOffset]
        : kAminoAcidCompositions['X' - kAAOffset];
}

} // namespace wasabi::aa