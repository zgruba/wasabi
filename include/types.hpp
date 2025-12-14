#pragma once
#include <string>
#include <sstream>
#include <vector>
#include <array>

namespace wasabi {

using comp_t = int16_t;

struct Composition {
    // Elements: C, H, N, O, S, P, Se
    std::array<comp_t, 7> elements{{0,0,0,0,0,0,0}};

    // ---- Constructors ----
    constexpr Composition() noexcept = default;

    constexpr explicit Composition(const std::array<comp_t,7>& elems) : elements(elems) {}

    constexpr Composition(const comp_t c, const comp_t h, const comp_t n, const comp_t o,
                         const comp_t s, const comp_t p, const comp_t se): elements{{c, h, n, o, s, p, se}} {}

    constexpr Composition(const Composition& other) noexcept {
        elements = other.elements;
    }

    Composition& operator=(const Composition& other) noexcept {
        if (this != &other)
            elements = other.elements;
        return *this;
    }

    Composition& operator+=(const Composition& other) noexcept {
        for (std::size_t i = 0; i < elements.size(); ++i)
            elements[i] += other.elements[i];
        return *this;
    }

    Composition& operator-=(const Composition& other) noexcept {
        for (std::size_t i = 0; i < elements.size(); ++i)
            elements[i] -= other.elements[i];
        return *this;
    }

    Composition& operator*=(int factor) noexcept {
        for (auto& e : elements)
            e *= factor;
        return *this;
    }

    [[nodiscard]] std::string to_string() const {
        static constexpr std::array<const char*,7> symbols{{"C","H","N","O","S","P","Se"}};
        std::ostringstream out;
        for (std::size_t i = 0; i < symbols.size(); ++i)
            if (elements[i] != 0) out << symbols[i] << elements[i];
        return out.str();
    }
};

inline Composition operator+(Composition a, const Composition& b) noexcept {
    return a += b;
}

inline Composition operator-(Composition a, const Composition& b) noexcept {
    return a -= b;
}

inline Composition operator*(Composition a, int f) noexcept {
    return a *= f;
}

using mz_t = double;
using rt_t = double;
using intensity_t = double;

struct PeakRow {
    double mz{};
    double rt{};        // seconds (OpenMS uses seconds by default)
    double intensity{};
};

struct FeatureRow {
    uint64_t feature_id{};
    double mz{};
    double min_mz{};
    double max_mz{};
    double rt{};
    double min_rt{};
    double max_rt{};
    double intensity{};
};

struct FastaRecord {
    std::string id;        // header line without '>' + seq {k} where k is k-th trypsin cut
    //std::string description; // remaining header text after first token
    std::string sequence;  // raw sequence (no whitespace)
    bool is_decoy;
};

using PeakTable = std::vector<PeakRow>;
using FeatureTable = std::vector<FeatureRow>;
using FastaRecords = std::vector<FastaRecord>;

inline std::ostream& operator<<(std::ostream& os, const PeakRow& row) {
    return os << "PeakRow(mz=" << row.mz
              << ", rt=" << row.rt
              << ", intensity=" << row.intensity << ")";
}

inline std::ostream& operator<<(std::ostream& os, const FeatureRow& row) {
    return os << "FeatureRow(id=" << row.feature_id
              << ", mz=" << row.mz
              << ", min_mz=" << row.min_mz
              << ", max_mz=" << row.max_mz
              << ", rt=" << row.rt
              << ", min_rt=" << row.min_rt
              << ", max_rt=" << row.max_rt
              << ", intensity=" << row.intensity << ")";
}

inline std::ostream& operator<<(std::ostream& os, const FastaRecord& rec) {
    return os << "FastaRecord(id=" << rec.id
             // << ", desc=" << rec.description
              << ", seq=" << rec.sequence << ")";
}

inline std::ostream& operator<<(std::ostream& os, const PeakTable& table) {
    os << "PeakTable[\n";
    for (const auto& row : table) {
        os << "  " << row << "\n";
    }
    return os << "]";
}

inline std::ostream& operator<<(std::ostream& os, const FeatureTable& table) {
    os << "FeatureTable[\n";
    for (const auto& row : table) {
        os << "  " << row << "\n";
    }
    return os << "]";
}

inline std::ostream& operator<<(std::ostream& os, const FastaRecords& records) {
    os << "FastaRecords[\n";
    for (const auto& rec : records) {
        os << "  " << rec << "\n";
    }
    return os << "]";
}

}
