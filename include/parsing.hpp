#pragma once

#include <string>
#include <mutex>
#include <thread>
#include <unordered_set>

#include "types.hpp"
#include "wasabi_db.hpp" // todo: get rid of that

namespace wasabi {

namespace fasta {


/// @brief Read a FASTA file into a list of protein records.
///
/// This function:
///   - Parses standard FASTA format (">" header + sequence lines)
///   - Extracts only the identifier (first token in header)
///   - Uppercases all amino acids
///   - Strips whitespace
///   - Does not apply any digestion or create decoys
///
/// @param file_path Path to the FASTA file.
/// @return Vector of FastaRecord objects, each containing {id, sequence, is_decoy=false}.
FastaRecords read_fasta(const std::string& file_path);

/// @brief Read a FASTA file and perform parallel trypsin digestion.
///
/// This routine:
///   - Loads all proteins from FASTA
///   - Digests each in parallel using trypsin rules (K/R not followed by P)
///   - Generates:
///         • Fully-tryptic peptides
///         • Optional semi-tryptic peptides (MaxQuant style)
///   - Enforces peptide length bounds (min_len–max_len)
///   - Supports up to `missed_cleavages` missed cuts
///   - Optionally generates MaxQuant-type decoy peptides
///     (reverse internal residues, termini unchanged)
///
/// @param fasta_path Path to FASTA file.
/// @param missed_cleavages Maximum allowed missed cleavages.
/// @param min_len Minimum peptide length.
/// @param max_len Maximum peptide length.
/// @param allow_semi Whether to generate semi-tryptic peptides.
/// @param allow_decoys If true, generate MaxQuant-style decoys.
/// @return Vector of digested peptide records (id, sequence, decoy flag).
FastaRecords read_fasta_with_trypsin(const std::string& fasta_path,
                                     int missed_cleavages,
                                     size_t min_len,
                                     size_t max_len,
                                     bool allow_semi,
                                     bool allow_decoys);


/// @brief Identify all trypsin cleavage sites in a protein sequence.
///
/// Trypsin rule implemented: (K|R)(!P)
///
/// Returned vector includes:
///   - 0 as the first cut site
///   - All valid internal cuts
///   - seq.size() as final boundary
///
/// @param seq Protein sequence.
/// @return Vector of indices marking cleavage boundaries.
std::vector<size_t> trypsin_cleavage_sites(const std::string& seq);

/// @brief Create a MaxQuant-style decoy by reversing internal residues.
///
/// For sequences >2 aa:
///   • Keep N-terminal and C-terminal residues fixed
///   • Reverse all internal residues
///
/// For very short sequences (<=2 aa), returns the input unchanged.
///
/// @param seq Original peptide or protein sequence.
/// @return Decoy sequence.
std::string make_maxquant_decoy(const std::string& seq);

/// @brief Digest a single protein into tryptic and semi-tryptic peptides.
///
/// This helper:
///   - Applies trypsin cuts with missed cleavages
///   - Generates fully and optionally semi-tryptic peptides
///   - Enforces min/max peptide length
///   - Ensures uniqueness (duplicates removed)
///   - Optionally generates matching decoys
///   - Appends results into `out` (thread-safe)
///
/// @param prot Input protein record.
/// @param out Output list where generated peptides are appended.
/// @param seen
/// @param missed_cleavages Allowed missed cleavages.
/// @param min_len Minimum peptide length.
/// @param max_len Maximum peptide length.
/// @param allow_semi Enable semi-tryptic digestion.
/// @param generate_decoys Enable MaxQuant-style decoys.
void digest_one_protein_dual(const FastaRecord& prot,
                                    std::vector<FastaRecord>& out,
                                    std::unordered_set<std::string>& seen,
                                    int missed_cleavages,
                                    size_t min_len,
                                    size_t max_len,
                                    bool allow_semi,
                                    bool generate_decoys);

} // namespace fasta

namespace mass {

/**
 * @brief Represents a feature's isotopic or peak distribution in m/z–intensity form.
 *
 * A FeatureSpec stores:
 *   - a textual identifier, and
 *   - a vector of m/z–intensity pairs describing the feature's spectral shape.
 *
 * Each entry in the distribution corresponds to a single peak, containing:
 *   - an m/z position, and
 *   - a normalized intensity or probability contribution.
 */
struct FeatureSpec {
    std::string id;   ///< Unique string identifier of the feature
    std::vector<db::DistributionElement<mz_t, intensity_t>>
        distribution; ///< Sorted m/z–intensity or m/z–probability distribution
};


/**
 * @brief Compact, array-based representation of a feature's isotopic distribution.
 *
 * CompactFeatureSpec stores the same conceptual information as FeatureSpec, but in a
 * memory-efficient form suitable for lookup, serialization,
 * and large-scale spectral databases.
 *
 * The structure contains:
 *   - a numeric identifier,
 *   - the number of isotopic peaks,
 *   - contiguous arrays of m/z values and their normalized intensities,
 *   - the intensity-weighted expected m/z of the distribution.
 *
 *
 * Example usage:
 *   - Storing theoretical isotope patterns in compressed database blocks
 *   - Rapid iteration during scoring
 */
struct CompactFeatureSpec {
    uint64_t id{};                     ///< Numeric feature identifier
    db::distribution_size_t size{};    ///< Number of peaks in the distribution
    std::unique_ptr<mz_t[]> values;    ///< Array of m/z values (length = size)
    std::unique_ptr<intensity_t[]> probabilities; ///< Normalized intensities (length = size)
    mz_t expected = 0.0;               ///< Intensity-weighted mean m/z of the distribution
    rt_t rt = 0.0;


    /// @brief Constructs an empty, uninitialized distribution.
    CompactFeatureSpec() = default;

    /**
     * @brief Constructs a feature distribution from existing arrays.
     *
     * @param id_       Numeric identifier of the feature
     * @param n         Number of m/z–intensity pairs
     * @param vals      Pointer to array of m/z values (ownership transferred)
     * @param probs     Pointer to array of normalized intensities (ownership transferred)
     * @param expected_ Precomputed expected m/z for the distribution
     * @param rt_       Retention time extracted from the feature file
     */
    CompactFeatureSpec(const uint64_t id_,
                 db::distribution_size_t n,
                 std::unique_ptr<mz_t[]> vals,
                 std::unique_ptr<intensity_t[]> probs,
                 const mz_t expected_,
                 const rt_t rt_) noexcept
        : id(id_),
          size(n),
          values(std::move(vals)),
          probabilities(std::move(probs)),
          expected(expected_),
          rt(rt_){}
};



/// @brief Load MS1 peaks from an mzML file.
///
/// Behavior:
///   - Loads an OpenMS::MSExperiment
///   - Extracts only MS1 spectra
///   - Flattens all MS1 peaks into a simple table:
///         { m/z, rt, intensity }
///
/// @param file_path Path to an mzML file.
/// @return PeakTable containing all MS1 peaks.
PeakTable load_mzml(const std::string& file_path);

/// @brief Load features from a featureXML file using OpenMS.
///
/// Extracts for each feature:
///   - Unique ID
///   - m/z, RT, intensity
///   - Convex hull bounds (min/max m/z, min/max RT)
///     if convex hulls exist
///
/// @param file_path Path to a featureXML file.
/// @return FeatureTable containing one row per OpenMS::Feature.
FeatureTable load_featurexml(const std::string& file_path);

/// @brief Generate a single feature's normalized isotopic-like distribution.
///
/// For a feature:
///   - Selects all MS1 points within its convex hull bounds
///   - Bins m/z values if requested
///   - Sums intensities within each bin
///   - Normalizes intensities to 1.0
///   - Computes expected m/z (intensity-weighted)
///   - Returns a FeatureSpec containing id + list of {mz, prob}
///
/// @param features FeatureTable containing extracted features.
/// @param points PeakTable containing MS1 peaks.
/// @param bin_mz If true, bin m/z values.
/// @param bin_width Width used for m/z binning.
/// @param threads Number of threads for parallel processing.
/// @return Vector of FeatureSpec objects, one per feature.
std::vector<FeatureSpec> generate_feature_entries(const FeatureTable& features,
                                                const PeakTable& points,
                                                bool bin_mz = false,
                                                double bin_width = 0.01,
                                                int threads = 10);

/// @brief Generate a single feature distribution in a compact form.
///
/// Same logic as generate_feature_entry(), except:
///   - Stores the distribution in raw arrays (unique_ptr)
///   - Stores the distribution size explicitly
///   - Stores expected m/z
///   - Uses numeric feature_id directly
///
/// @param features FeatureTable containing extracted features.
/// @param points PeakTable containing MS1 peaks.
/// @param bin_mz Enable m/z binning.
/// @param bin_width Bin width if binning enabled.
/// @param threads Number of threads.
/// @return Vector of CompactFeatureSpec objects, sorted by expected m/z.
std::vector<CompactFeatureSpec> generate_compact_feature_entries(
    const FeatureTable& features,
    const PeakTable& points,
    bool bin_mz,
    double bin_width,
    int threads);

} // namespace mass
} // namespace wasabi
