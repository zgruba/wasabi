#include "parsing.hpp"

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <stdexcept>
#include <fstream>
#include <cctype>
#include <mutex>
#include <vector>
#include <string>
#include <algorithm>
#include <future>

#include <omp.h>

#include "wasabi_db.hpp"
#include "utils.hpp"

namespace wasabi {
// ----------------- Utilities -----------------
namespace util {
    inline char to_upper(const char c) {
        return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }

    inline bool is_space(const char c) {
        return std::isspace(static_cast<unsigned char>(c));
    }
}

// ----------------- FASTA -----------------
namespace fasta {

bool is_valid(const char c) {
    switch (util::to_upper(c)) {
        case 'A': case 'R': case 'N': case 'D': case 'C':
        case 'E': case 'Q': case 'G': case 'H': case 'I':
        case 'L': case 'K': case 'M': case 'F': case 'P':
        case 'S': case 'T': case 'W': case 'Y': case 'V':
        case 'B': case 'J': case 'O': case 'U': case 'X': case 'Z':
        case '*': case '-':
            return true;
        default:
            return false;
    }
}
// --------------------------------------------------------------------
// FASTA reader (assumes valid FASTA)
// --------------------------------------------------------------------
FastaRecords read_fasta(const std::string& file_path) {
    std::ifstream in(file_path);
    if (!in)
        throw std::runtime_error("Cannot open FASTA: " + file_path);

    FastaRecords records;
    FastaRecord cur;
    std::string line;

    auto flush = [&]() {
        if (!cur.id.empty() && !cur.sequence.empty()) {
            records.push_back(std::move(cur));
            cur = FastaRecord{};
        }
    };

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            flush();
            std::string header = line.substr(1);
            if (auto pos = header.find_first_of(" \t"); pos == std::string::npos) {
                cur.id = header;
            } else {
                cur.id = header.substr(0, pos);
            }
            cur.is_decoy = false;
        } else {
            for (char c : line) {
                if (util::is_space(c)) continue;
                cur.sequence.push_back(util::to_upper(c));
            }
        }
    }
    flush();
    return records;
}

// --------------------------------------------------------------------
// Trypsin cleavage sites: (K|R)(?!P)
// --------------------------------------------------------------------
inline std::vector<size_t> trypsin_cleavage_sites(const std::string& seq) {
    std::vector<size_t> cuts;
    //cuts.reserve(seq.size() / 10 + 4);
    cuts.push_back(0);

    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        const char c = seq[i];
        const char n = seq[i + 1];
        if ((c == 'K' || c == 'R') && n != 'P')
            cuts.push_back(i + 1);
    }

    cuts.push_back(seq.size());
    return cuts;
}

// --------------------------------------------------------------------
// MaxQuant decoy: reverse internal residues (keep termini fixed)
// --------------------------------------------------------------------
inline std::string make_maxquant_decoy(const std::string& seq) {
    if (seq.size() <= 2) return seq;
    std::string rev = seq;
    std::reverse(rev.begin() + 1, rev.end() - 1);
    return rev;
}

inline void digest_one_protein_dual(const FastaRecord& prot,
                                    std::vector<FastaRecord>& out,
                                    std::unordered_set<std::string>& seen,
                                    const int missed_cleavages,
                                    const size_t min_len,
                                    const size_t max_len,
                                    const bool allow_semi,
                                    const bool generate_decoys)
{
    std::vector<FastaRecord> local;
    if (!generate_decoys) {
        const auto cuts = trypsin_cleavage_sites(prot.sequence);
        size_t k = 0;
        const std::string& seq = prot.sequence;

        for (size_t i = 0; i + 1 < cuts.size(); ++i) {
            for (int mc = 0; mc <= missed_cleavages; ++mc) {
                const size_t j = i + 1 + mc;
                if (j >= cuts.size()) break;

                const size_t tryptic_start = cuts[i];
                const size_t tryptic_end   = cuts[j];
                const size_t tryptic_len   = tryptic_end - tryptic_start;

                // ---------- Fully tryptic peptide ----------
                if (tryptic_len >= min_len && tryptic_len <= max_len) {
                    std::string pep = seq.substr(tryptic_start, tryptic_len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "{" + std::to_string(k++) + "}", pep, false });
                    }
                }

                if (!allow_semi) continue;
                // ---------- Semi-tryptic peptides (MaxQuant style) ----------
                // N-semi: internal K/R -> tryptic_end
                for (size_t start = tryptic_start + 1; start < tryptic_end; ++start) {
                    const size_t len = tryptic_end - start;
                    if (len < min_len || max_len < len) continue;

                    std::string pep = seq.substr(start, len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "{" + std::to_string(k++) + "}", pep, false });
                    }
                }

                // C-semi: tryptic_start -> internal K/R
                for (size_t end = tryptic_start + 1; end < tryptic_end; ++end) {
                    const size_t len = end - tryptic_start;
                    if (len < min_len || len > max_len) continue;

                    std::string pep = seq.substr(tryptic_start, len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "{" + std::to_string(k++) + "}", pep, false });
                    }
                }
            }
        }
    }
    // If generate_decoys == true apply similar logic for decoys
    else {
        std::string decoy_seq;
        decoy_seq = make_maxquant_decoy(prot.sequence);
        const auto decoy_cuts = trypsin_cleavage_sites(decoy_seq);
        size_t k = 0;
        for (size_t i = 0; i + 1 < decoy_cuts.size(); ++i) {
            for (int mc = 0; mc <= missed_cleavages; ++mc) {
                const size_t j = i + 1 + mc;
                if (j >= decoy_cuts.size()) break;

                const size_t tryptic_start = decoy_cuts[i];
                const size_t tryptic_end   = decoy_cuts[j];
                const size_t tryptic_len   = tryptic_end - tryptic_start;

                // ---------- Fully tryptic peptide ----------
                if (tryptic_len >= min_len && tryptic_len <= max_len) {
                    std::string pep = decoy_seq.substr(tryptic_start, tryptic_len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "_DECOY{" + std::to_string(k++) + "}", pep, true });
                    }
                }

                if (!allow_semi) continue;
                // ---------- Semi-tryptic peptides (MaxQuant style) ----------
                // N-semi: internal K/R -> tryptic_end
                for (size_t start = tryptic_start + 1; start < tryptic_end; ++start) {
                    const size_t len = tryptic_end - start;
                    if (len < min_len || max_len < len) continue;

                    std::string pep = decoy_seq.substr(tryptic_start, tryptic_len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "_DECOY{" + std::to_string(k++) + "}", pep, true });
                    }
                }

                // C-semi: tryptic_start -> internal K/R
                for (size_t end = tryptic_start + 1; end < tryptic_end; ++end) {
                    const size_t len = end - tryptic_start;
                    if (len < min_len || len > max_len) continue;

                    std::string pep = decoy_seq.substr(tryptic_start, tryptic_len);
                    if (seen.insert(pep).second) {
                        local.push_back({ prot.id + "_DECOY{" + std::to_string(k++) + "}", pep, true });
                    }
                }
            }
        }
    }

    // Thread-safe append
    static std::mutex out_mtx;
    {
        std::lock_guard<std::mutex> lock(out_mtx);
        out.insert(out.end(),
                   std::make_move_iterator(local.begin()),
                   std::make_move_iterator(local.end()));
    }
}



// --------------------------------------------------------------------
// Master orchestrator (with decoys controlled by flag)
// --------------------------------------------------------------------
FastaRecords read_fasta_with_trypsin(const std::string& fasta_path,
                                     const int missed_cleavages,
                                     const size_t min_len,
                                     const size_t max_len,
                                     const bool allow_semi,
                                     const bool allow_decoys)
{
    const auto proteins = read_fasta(fasta_path);
    FastaRecords out;
    std::unordered_set<std::string> seen;

    // Parallelizing these loops gave very little benefit since here wer can cary over 'seen' between loops more easily.
    for (auto& prot : proteins) {
        digest_one_protein_dual(prot, out, seen, missed_cleavages, min_len, max_len, allow_semi, false);
    }
    if (allow_decoys) {
        for (auto& prot : proteins) {
            digest_one_protein_dual(prot, out, seen, missed_cleavages, min_len, max_len, allow_semi, true);
        }
    }

    out.shrink_to_fit();
    return out;
}

}


// ----------------- OpenMS -----------------
namespace mass {
PeakTable load_mzml(const std::string& file_path) {
    OpenMS::MSExperiment exp;
    try {
        OpenMS::MzMLFile mzml;
        mzml.load(file_path, exp);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string{"Failed to load mzML '"} + file_path + "': " + e.what());
    }

    PeakTable rows;
    for (const OpenMS::MSSpectrum& spectrum : exp) {
        if (spectrum.getMSLevel() != 1) continue; // MS1 only
        const double rt = spectrum.getRT();
        for (const OpenMS::Peak1D& p : spectrum) {
            rows.push_back({p.getMZ(), rt, p.getIntensity()});
        }
    }

    return rows;
}

FeatureTable load_featurexml(const std::string& file_path) {
    OpenMS::FeatureMap fmap;
    try {
        OpenMS::FeatureXMLFile fxml;
        fxml.load(file_path, fmap);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string{"Failed to load featureXML '"} + file_path + "': " + e.what());
    }

    FeatureTable rows;
    rows.reserve(fmap.size());

    for (const OpenMS::Feature& f : fmap) {
        const auto uid = static_cast<uint64_t>(f.getUniqueId());
        const double mz = f.getMZ();
        const double rt = f.getRT();
        const double intensity = f.getIntensity();

        double min_mz = std::numeric_limits<double>::max();
        double max_mz = std::numeric_limits<double>::min();
        double min_rt = std::numeric_limits<double>::max();
        double max_rt = std::numeric_limits<double>::min();

        // Incorporate convex hulls if they exist
        if (const auto& hulls = f.getConvexHulls(); !hulls.empty()) {
            for (const auto& hull : hulls) {
                for (const auto& point : hull.getHullPoints()) {
                    // point[0] - rt, point[1] - mz
                    min_rt = std::min(min_rt, point[0]);
                    min_mz = std::min(min_mz, point[1]);
                    max_rt = std::max(max_rt, point[0]);
                    max_mz = std::max(max_mz, point[1]);
                }
            }
        }
        // Push row with convex hull bounds
        rows.push_back({uid, mz, min_mz, max_mz, rt, min_rt, max_rt, intensity});
    }

    return rows;
}

// ---------------- Generate Entries for each feature ----------------
FeatureSpec generate_feature_entry(const FeatureRow& feat, const PeakTable& points,
                                                const bool bin_mz, const double bin_width) {
    // 1. Select points inside the feature convex hull bounds
    std::vector<PeakRow> subset;
    for (const auto& p : points) {
        if (p.mz >= feat.min_mz && p.mz <= feat.max_mz && p.rt >= feat.min_rt && p.rt <= feat.max_rt) {
            subset.push_back(p);
        }
    }
    if (subset.empty()) return {};

    // 2. Group by m/z or bin
    std::map<double,  double> grouped; // m/z -> summed intensity
    for (const auto& p : subset) {
        double key = p.mz;
        if (bin_mz) {
            key = std::floor(p.mz / bin_width) * bin_width;
        }
        grouped[key] += (p.intensity);
    }

    // 3. Normalize intensities
    double total_intensity = 0.0;
    for (const auto& [mz, intensity] : grouped) {
        total_intensity += intensity;
    }
    if (total_intensity > 0.0) {
        for (auto& [mz, intensity] : grouped) {
            intensity /= total_intensity;
        }
    }


    // 4. Compute expected m/z (weighted average)
    double expected_mz = 0.0;
    for (const auto& [mz, intensity] : grouped) {
        expected_mz += mz * intensity;
    }


    // 5. Create Entry
    FeatureSpec e;
    e.distribution.clear();
    e.id = std::to_string(feat.feature_id);

    for (const auto& [mz, intensity] : grouped) {
        e.distribution.push_back({mz, static_cast<double>(intensity)});
    }

    return e;
}

std::vector<FeatureSpec> generate_feature_entries(
    const FeatureTable& features,
    const PeakTable& points,
    const bool bin_mz,
    const double bin_width,
    const int threads)
{
    // Each thread stores its own results
    std::vector<std::vector<FeatureSpec>> local_results(static_cast<size_t>(std::max(1, threads)));

#pragma omp parallel num_threads(threads) default(none) shared(local_results, features, points, bin_mz, bin_width)
    {
        const int tid = omp_get_thread_num();
        auto& local = local_results[static_cast<size_t>(tid)];

#pragma omp for schedule(static)
        for (const auto& feature : features) {
            local.push_back(generate_feature_entry(feature, points, bin_mz, bin_width));
        }
    }

    // Merge all thread results into one vector
    std::vector<FeatureSpec> results;
    size_t total_size = 0;
    for (const auto& v : local_results) {
        total_size += v.size();
    }
    results.reserve(total_size);

    for (auto& v : local_results) {
        std::ranges::move(v, std::back_inserter(results));
    }

    return results;
}

// ---------------- Generate Entries (compact) for each feature ----------------
CompactFeatureSpec generate_compact_feature_entry(
    const FeatureRow& feat,
    const PeakTable& points,
    const bool bin_mz,
    const double bin_width)
{
    // 1. Select relevant peaks within feature bounds
    std::vector<PeakRow> subset;
    subset.reserve(points.size());
    for (const auto& p : points) {
        if (p.mz >= feat.min_mz && p.mz <= feat.max_mz &&
            p.rt >= feat.min_rt && p.rt <= feat.max_rt) {
            subset.push_back(p);
        }
    }

    if (subset.empty()) {
        return {}; // empty CompactFeatureSpec
    }

    // 2. Group by m/z (bin if requested)
    std::map<double, double> grouped;
    for (const auto& p : subset) {
        double key = bin_mz
            ? std::floor(p.mz / bin_width) * bin_width
            : p.mz;
        grouped[key] += p.intensity;
    }

    if (grouped.empty()) {
        return {};
    }

    // 3. Normalize intensities
    double total_intensity = 0.0;
    for (const auto& [_, intensity] : grouped) {
        total_intensity += intensity;
    }

    if (total_intensity <= 0.0) {
        return {};
    }

    for (auto& [_, intensity] : grouped) {
        intensity /= total_intensity;
    }

    // 4. Compute expected m/z (intensity-weighted average)
    double weighted_sum = 0.0;
    for (const auto& [mz, intensity] : grouped) {
        weighted_sum += mz * intensity;
    }

    const mz_t expected = static_cast<mz_t>(weighted_sum);

    // 5. Allocate arrays for values and probabilities
    const auto n = static_cast<db::distribution_size_t>(grouped.size());
    auto values = std::make_unique<mz_t[]>(n);
    auto probs = std::make_unique<intensity_t[]>(n);

    size_t i = 0;
    for (const auto& [mz, intensity] : grouped) {
        values[i] = static_cast<mz_t>(mz);
        probs[i]  = static_cast<intensity_t>(intensity);
        ++i;
    }

    // 6. Calculate mean rt
    const auto rt = feat.rt;

    // 7. Construct and return CompactFeatureSpec
    return CompactFeatureSpec{
        feat.feature_id,
        n,
        std::move(values),
        std::move(probs),
        expected,
        rt
    };
}


std::vector<CompactFeatureSpec> generate_compact_feature_entries(
    const FeatureTable& features,
    const PeakTable& points,
    const bool bin_mz,
    const double bin_width,
    const int threads)
{
    if (features.empty()) {
        return {};
    }

    const int num_threads = std::max(1, threads);
    std::vector<std::vector<CompactFeatureSpec>> local_results(static_cast<size_t>(num_threads));

#pragma omp parallel num_threads(num_threads) default(none) shared(local_results, features, points, bin_mz, bin_width, threads)
    {
        const int tid = omp_get_thread_num();
        auto& local = local_results[static_cast<size_t>(tid)];

#pragma omp for schedule(static)
        for (const auto & feat : features) {
            if (auto spec = generate_compact_feature_entry(feat, points, bin_mz, bin_width); spec.size > 0) {
                local.push_back(std::move(spec));
            }
        }
    }

    // Combine all thread-local vectors efficiently
    size_t total_size = 0;
    for (const auto& v : local_results) {
        total_size += v.size();
    }

    std::vector<CompactFeatureSpec> results;
    results.reserve(total_size);

    for (auto& v : local_results) {
        std::ranges::move(v, std::back_inserter(results));
    }

    // Sort by expected m/z for deterministic downstream analysis
    std::ranges::sort(results,
                      [](const CompactFeatureSpec& a, const CompactFeatureSpec& b) {
                          return a.expected < b.expected;
                      });

    return results;
}

}
}
