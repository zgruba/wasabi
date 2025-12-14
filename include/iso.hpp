#pragma once
#include <string>
#include <vector>
#include <map>
#include <memory>

#include <IsoSpec++/isoSpec++.h>
#include <boost/intrusive_ptr.hpp>

#include "types.hpp"
#include "utils.hpp"
#include "ptm.hpp"

namespace wasabi::iso {
struct ModifiedPeptide;
struct AnnotationInfo;

using entry_id_t = uint64_t;
using seq_id_t = uint32_t;
using ptm_id_t = ptm::common_ptms_size_t;
using ion_id_t = uint8_t;
using charge_t = int8_t;
using ptm_count_t = uint16_t;
using dist_size_t = uint16_t;
using multi_t = int8_t;

/**
 * @brief Describes ionisation properties used for m/z computation.
 *
 * An ionisation mode defines how a neutral peptide becomes an ion.
 * It specifies:
 *   - the resulting charge state (e.g. +1 for protonation, +2 for doubly charged),
 *   - the mass shift contributed by the adduct or ionising species.
 *
 * Examples:
 *   - Protonation: charge = +1, mass_shift = +1.007276 (mass of H⁺)
 *   - Sodium adduct: charge = +1, mass_shift = +22.989218 (mass of Na⁺)
 *
 * This information is used when converting neutral peptide masses
 * into the final m/z value used in isotopic distribution calculations.
 */
struct IonisationMode {
    charge_t charge;        ///< Ion charge state (positive or negative)
    double mass_shift;      ///< Mass added by the ionising species before dividing by charge
};


/**
 * @brief A fully specified peptide instance with all modifications, adducts,
 *        and charge state, used during isotopic distribution generation.
 *
 * A PeptideEntry represents a unique “chemical species” produced during
 * peptide enumeration. It binds together:
 *   - a peptide sequence (via ModifiedPeptide) including PTMs,
 *   - an ionisation/adduct type,
 *   - the resulting charge state,
 *   - the full elemental composition after applying all PTMs and the adduct,
 *   - decoy/target labelling.
 *   - predicted rt (default NaN)
 *
 * PeptideEntry objects are stored in sorted collections, so comparison
 * operators order entries primarily by elemental composition,
 * then by charge, then by adduct name.
 *
 * These entries serve as the input for computing isotopic envelopes.
 */
struct PeptideEntry {
    const std::string* id;                             ///< Peptide string ID (sequence or external ID)
    boost::intrusive_ptr<ModifiedPeptide> peptide;      ///< Pointer to sequence + PTMs
    charge_t charge;                                    ///< Final ion charge after ionisation
    const std::string* adduct_name;                     ///< Human-readable adduct/ionisation label
    Composition total_composition;                      ///< Full elemental composition of the ionised species
    bool is_decoy = false;                              ///< Whether this peptide is a decoy entry
    double rt = std::numeric_limits<double>::quiet_NaN();

    PeptideEntry(const std::string* id_,
                 const boost::intrusive_ptr<ModifiedPeptide>& pep_,
                 charge_t c_,
                 const std::string* adduct_,
                 const Composition& comp_,
                 bool decoy)
       : id(id_), peptide(pep_), charge(c_),
         adduct_name(adduct_), total_composition(comp_), is_decoy(decoy)
    {}

    bool operator<(const PeptideEntry& other) const {
        if (total_composition.elements < other.total_composition.elements)
            return true;
        if (other.total_composition.elements < total_composition.elements)
            return false;

        if (charge < other.charge)
            return true;
        if (other.charge < charge)
            return false;

        if (*adduct_name < *other.adduct_name)
            return true;
        if (*adduct_name > *other.adduct_name)
            return false;

        return static_cast<int>(is_decoy) < static_cast<int>(other.is_decoy);
    }

    bool operator==(const PeptideEntry& other) const {
        return total_composition.elements == other.total_composition.elements &&
               charge == other.charge &&
               *adduct_name == *other.adduct_name;
    }
};


/**
 * @brief Describes the isotopic distribution metadata for a peptide entry.
 *
 * Each IsotopicInfo record corresponds to a single PeptideEntry after
 * ionisation and composition evaluation. It stores:
 *   - the peptide ID,
 *   - the expected monoisotopic mass or m/z value,
 *   - the number of isotopic peaks to store,
 *   - the adduct name (e.g., "protonation", "sodium"),
 *   - the final charge state after all ionisation effects.
 *
 * The structure is designed as a lightweight header used before reading
 * the full isotopic distribution from disk or memory.
 */
struct IsotopicInfo {
    entry_id_t id;                   ///< Peptide index linking back to PeptideEntry
    mz_t expected;                   ///< Expected monoisotopic neutral mass or m/z (depending on context)
    dist_size_t distribution_size;   ///< Number of isotopic peaks in the computed envelope
    const std::string* adduct_name;  ///< Ionisation/adduct type
    charge_t total_charge;           ///< Final total charge of the ion
};


/**
 * @brief A compact on-disk header describing where an isotopic distribution is stored.
 *
 * HeaderRecord is written into the distribution index file and allows
 * fast random access into the isotopic distribution data block.
 *
 * It stores:
 *   - the peptide ID,
 *   - the expected monoisotopic mass or m/z,
 *   - the number of peaks in the distribution,
 *   - the byte offset in the associated data file.
 */
struct HeaderRecord {
    entry_id_t id = 0;           ///< Peptide ID
    mz_t expected = 0.0;         ///< Monoisotopic expected mass/m/z
    dist_size_t size = 0;        ///< Size of the isotopic envelope
    uint64_t offset = 0;         ///< Byte offset into distribution storage

    HeaderRecord() = default;

    HeaderRecord(entry_id_t id_, mz_t expected_, dist_size_t size_, uint64_t offset_)
        : id(id_), expected(expected_), size(size_), offset(offset_) {}
};


/**
 * @brief Computes total elemental formula by adding modification formulas to a peptide base formula.
 *
 * @param peptide_formula Base chemical formula of peptide (e.g. "C15H21N4O11S1P1").
 * @param mods Map of modification names to their applied count.
 * @return Combined total formula string.
 *
 * @note Formula parsing is handled element-wise, not by string concatenation.
 */
std::string compute_total_formula(const std::string& peptide_formula,
                                  const std::map<std::string, uint16_t>& mods);

/**
 * @brief Computes a single m/z value from a neutral mass given ionisation mode and charge multiplier.
 *
 * @param mass Neutral mass (Da)
 * @param ion_info Ionisation mode (charge and shift)
 * @param multiplier Charge multiplier
 * @return Computed m/z value
 *
 * @note Uses formula: (mass + shift * multiplier) / |charge * multiplier|
 */
double compute_mz(double mass, const IonisationMode& ion_info, int multiplier);

/**
 * @brief Generates isotopic envelope (masses and probabilities) for a given chemical formula.
 *
 * @param pep Peptide sequence
 * @param threshold Total probability cutoff for isotope inclusion (e.g. 0.999)
 * @param[out] mzs Output vector of corresponding values
 * @param[out] probs Output vector of corresponding probabilities
 *
 * @details Uses IsoSpec++’s IsoThresholdGenerator for isotope generation.
 */
void generate_isotopic_distribution(const PeptideEntry& pep, double threshold,
                                    std::vector<double>& mzs, std::vector<double>& probs);

/**
 * @brief Represents a peptide variant with applied chemical modifications.
 *
 * This structure encapsulates a peptide sequence together with all post-translational
 * and chemical modifications applied to it. It acts as the canonical “modified peptide”
 * object used throughout enumeration, mass calculation, and isotopic distribution steps.
 *
 * A ModifiedPeptide stores:
 *   - a non-owning pointer to the original amino-acid sequence (from FASTA),
 *   - residue-level modifications (e.g., Oxidation on M, Carbamidomethyl on C),
 *   - terminal modifications on the N-terminus or C-terminus,
 *   - an intrusive reference counter for efficient deduplication and sharing.
 *
 * Residue modifications are stored in a compact array `residue_mods_ptr` of pairs:
 * `(count, ptm_index)` where:
 *   - `count` is how many times the modification occurs at that residue
 *     (typically 0 or 1, but higher counts allow stacked modifications),
 *   - `ptm_index` is an index into the global `common_ptms_k` table.
 *
 * The array is laid out in order of occurrence along the sequence, and its size is
 * stored in `residue_mods_size`. If no residues are modified, this size is zero.
 *
 * ### Terminal modifications
 * The fields `n_term_mod_idx` and `c_term_mod_idx` hold the PTM indices for
 * N-terminal and C-terminal modifications. If a terminus is unmodified, the value is
 * `ptm::no_op_ptm_idx`. Convenience methods `has_n_term()` and `has_c_term()`
 * expose these flags, and getters return the PTM descriptors.
 *
 * ### Memory and ownership
 * The object owns its allocated array of residue modifications and cleans it up
 * in the destructor. It is managed by `boost::intrusive_ptr`, so copying is explicit
 * and reference counting is performed through `intrusive_ptr_add_ref` /
 * `intrusive_ptr_release`.
 *
 * This design allows ModifiedPeptide instances to be cheaply shared between many
 * PeptideEntry objects while preventing accidental duplication of heavy modification
 * data structures.
 */
struct ModifiedPeptide {
    const std::string* sequence;                                  // non-owning view into FASTA record
    std::pair<uint8_t, ptm::common_ptms_size_t>* residue_mods_ptr; // (number of applied modifications, PTM idx)
    ptm::common_ptms_size_t residue_mods_size;
    ptm::common_ptms_size_t c_term_mod_idx = ptm::no_op_ptm_idx;
    ptm::common_ptms_size_t n_term_mod_idx = ptm::no_op_ptm_idx;

    ModifiedPeptide(const std::string* sequence_,
                    std::vector<std::pair<uint8_t, ptm::common_ptms_size_t>>&& residueMods,
                    const ptm::common_ptms_size_t c_term_mod_idx_, const ptm::common_ptms_size_t n_term_mod_idx_) :
        sequence(sequence_),  c_term_mod_idx(c_term_mod_idx_), n_term_mod_idx(n_term_mod_idx_), intrusive_counter_(0) {
        residue_mods_size = static_cast<ptm::common_ptms_size_t>(residueMods.size());
        residue_mods_ptr = new std::pair<uint8_t, ptm::common_ptms_size_t>[residue_mods_size];
        for (size_t i = 0; i < residue_mods_size; i++) {
            residue_mods_ptr[i] = residueMods[i];
        }
    }

    ~ModifiedPeptide() {
        delete[] residue_mods_ptr;
    };

    ModifiedPeptide(const ModifiedPeptide& other) {
        sequence = other.sequence;
        residue_mods_size = other.residue_mods_size;
        residue_mods_ptr = new std::pair<uint8_t, ptm::common_ptms_size_t>[residue_mods_size];
        for (size_t i = 0; i < residue_mods_size; i++) {
            residue_mods_ptr[i] = other.residue_mods_ptr[i];
        }
        c_term_mod_idx = other.c_term_mod_idx;
        n_term_mod_idx = other.n_term_mod_idx;
        intrusive_counter_ = 0;
    }

    ModifiedPeptide(ModifiedPeptide&& other) noexcept {
        sequence = other.sequence;
        residue_mods_size = other.residue_mods_size;
        residue_mods_ptr = new std::pair<uint8_t, ptm::common_ptms_size_t>[residue_mods_size];
        for (size_t i = 0; i < residue_mods_size; i++) {
            residue_mods_ptr[i] = other.residue_mods_ptr[i];
        }
        c_term_mod_idx = other.c_term_mod_idx;
        n_term_mod_idx = other.n_term_mod_idx;
        intrusive_counter_ = other.intrusive_counter_;
    }

    ModifiedPeptide& operator=(const ModifiedPeptide& other) {
        if (this == &other) {
            return *this;
        }
        sequence = other.sequence;
        residue_mods_size = other.residue_mods_size;
        residue_mods_ptr = new std::pair<uint8_t, ptm::common_ptms_size_t>[residue_mods_size];
        for (size_t i = 0; i < residue_mods_size; i++) {
            residue_mods_ptr[i] = other.residue_mods_ptr[i];
        }
        c_term_mod_idx = other.c_term_mod_idx;
        n_term_mod_idx = other.n_term_mod_idx;
        intrusive_counter_ = 0;
        return *this;
    }

    ModifiedPeptide& operator=(ModifiedPeptide&& other) noexcept {
        sequence = other.sequence;
        residue_mods_size = other.residue_mods_size;
        residue_mods_ptr = new std::pair<uint8_t, ptm::common_ptms_size_t>[residue_mods_size];
        for (size_t i = 0; i < residue_mods_size; i++) {
            residue_mods_ptr[i] = other.residue_mods_ptr[i];
        }
        c_term_mod_idx = other.c_term_mod_idx;
        n_term_mod_idx = other.n_term_mod_idx;
        intrusive_counter_ = other.intrusive_counter_;
        return *this;
    }

    [[nodiscard]] bool has_n_term() const {return n_term_mod_idx != ptm::no_op_ptm_idx;}
    [[nodiscard]] bool has_c_term() const {return c_term_mod_idx != ptm::no_op_ptm_idx;}

    [[nodiscard]] const ptm::PTM& get_n_term() const { return ptm::common_ptms_k[n_term_mod_idx]; }
    [[nodiscard]] const ptm::PTM& get_c_term() const { return ptm::common_ptms_k[c_term_mod_idx]; }

    [[nodiscard]] ptm::common_ptms_size_t terminal_mods_count() const { return has_c_term() + has_n_term(); }

    friend void intrusive_ptr_add_ref(ModifiedPeptide* p) { ++p->intrusive_counter_; }
    friend void intrusive_ptr_release(ModifiedPeptide* p) {
        if (--p->intrusive_counter_ == 0) delete p;
    }

    bool operator==(const ModifiedPeptide& other) const {
        if (residue_mods_size == other.residue_mods_size
            && c_term_mod_idx == other.c_term_mod_idx
            && n_term_mod_idx == other.n_term_mod_idx
            && *sequence == *other.sequence) {
            for (size_t i = 1; i < residue_mods_size; i++) {
                if (residue_mods_ptr[i] != other.residue_mods_ptr[i]) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
private:
    uint32_t intrusive_counter_;
};

/**
 * @brief Generates all modified peptide variants (fixed and variable) from a peptide sequence.
 *
 * @param sequence Peptide sequence (e.g. "ACDEFGHIK")
 * @return List of ModifiedPeptide objects containing formulas and applied modifications.
 */
std::vector<ModifiedPeptide> generate_modified_peptides(std::string* sequence);

/**
 * @brief Generate all charge/adduct PeptideEntry variants for a single FASTA sequence,
 *        including all variable and fixed PTM combinations.
 *
 * This function:
 *   - Enumerates all ModifiedPeptide variants using generate_modified_peptides_site().
 *   - Computes chemical composition.
 *   - Creates PeptideEntry objects for every ionisation mode and charge multiplier.
 *   - Does NOT compute isotopic distributions.
 *
 * @param records FASTA-derived peptide entries
 * @param multipliers Charge state multipliers (e.g. {1, 2, 3})
 * @param threads number of threads to run the generation on
 * @param max_variable_mods
 * @return Vector of database entries representing simulated spectra.
 */

std::vector<PeptideEntry> generate_entries_from_fasta_records(
    const std::vector<FastaRecord>& records,
    const std::vector<multi_t>& multipliers,
    size_t max_variable_mods,
    int threads = 10);

/**
 * @brief Generates a single modified spectrum entry for a given peptide and modification map.
 *
 * @param fasta_id Peptide identifier from FASTA
 * @param fasta_sequence Peptide sequence
 * @param mods Applied modification counts
 * @param adduct_name Adduct ionisation mode (default: "protonation")
 * @param charge_multiplier Charge state multiplier (default: +1)
 * @return Corresponing PeptideEntry object.
 */
PeptideEntry generate_modified_entry(const std::string& fasta_id,
                                  const std::string& fasta_sequence,
                                  const std::map<std::string, uint16_t>& mods,
                                  const std::string& adduct_name,
                                  multi_t charge_multiplier);

/**
 * @brief Writes peptide entries and related metadata to a binary file.
 *
 * This binary format stores:
 *   - all peptide identifiers,
 *   - ModifiedPeptide descriptors,
 *   - total elemental compositions,
 *   - adducts and charges,
 *
 * The binary layout is optimized for fast loading and reproducible
 * annotation. Strings are deduplicated and written once.
 *
 * @param filename Path to output binary file.
 * @param fasta_records Original FASTA records (used for mapping).
 * @param entries Vector of PeptideEntry objects.
 */
void write_binary_info(const std::string& filename,
                       const FastaRecords& fasta_records,
                       const std::vector<PeptideEntry>& entries);

/**
 * @brief Read annotations written by write_binary_info() and reconstruct
 *        AnnotationInfo objects with shared string resources.
 *
 * Recreates:
 *   - sequence table
 *   - PTM name table
 *   - ion table
 *   - ModifiedPeptide including residue and terminal PTMs
 *   - chemical composition (recomputed)
 *
 * @param filename The binary file previously written with write_binary_info().
 * @return Vector of annotation records restored from disk.
 */
std::vector<AnnotationInfo> read_from_binary_info(const std::string& filename);


/**
* @brief Generate only the *summary* isotopic info (expected m/z) for sorting.
 *
 * Computes:
 *   - isotopic distribution until cumulative probability >= threshold
 *   - normalized probabilities
 *   - expected mass
 *   - expected m/z using adduct and charge
 *
 *
 * @param pep Peptide entry containing sequence, composition, adduct, charge.
 * @param threshold Relative intensity cutoff (0–1).
 * @param id Numerical index for cross-referencing.
 * @return IsotopicInfo containing peaks, expected value, and indexing.
 */
IsotopicInfo generate_isotopic_info(const PeptideEntry& pep, double threshold, size_t id);

/**
 * @brief Sorts peptides by their expected isotopic mass (or mz).
 *
 * Expected isotopic mass is computed using generate_isotopic_info() for
 * each peptide. Sorting is stable and ascending.
 *
 * @param peptides Input peptides.
 * @param threshold Intensity cutoff passed to isotopic calculation.
 * @return Sorted vector of IsotopicInfo objects.
 */
std::vector<IsotopicInfo>
sort_peptides_by_expected_value(const std::vector<PeptideEntry>& peptides,
                                double threshold);

/**
* @brief Write full isotopic distributions to a binary file.
 *
 * Layout:
 *   uint32 N
 *   HeaderRecord headers[N]
 *       each contains: id, expected_mz, distribution_size, offset
 *
 *   For each entry:
 *       mz_t    m/z values [distribution_size]
 *       intensity_t probabilities [distribution_size]
 *
 * @param filename Output file.
 * @param entries Source peptide entries.
 * @param isotopic_infos Precomputed isotopic information.
 * @param threshold Filtering threshold used during generation.
 */
void write_binary_distribution(const std::string& filename,
                               const std::vector<PeptideEntry>& entries,
                               const std::vector<IsotopicInfo>& isotopic_infos,
                               double threshold);

/**
 * @brief Reads full isotopic distributions from binary.
 *
 * Mostly intended for debugging and regression testing. Each record contains:
 *   - peptide index,
 *   - center-of-mass,
 *   - vector of m/z values,
 *   - vector of intensities.
 *
 * @param filename Binary file.
 * @return Vector of (index, expected_mass, values[], intensities[]) tuples.
 */
std::vector<std::tuple<size_t, mz_t, std::vector<mz_t>, std::vector<intensity_t>>>
read_binary_distribution(const std::string& filename); // mainly for testing

/**
 * @brief Shared string-storage for annotation objects.
 *
 * AnnotationInfo stores raw pointers into these vectors in order to
 * avoid string duplication. The shared_ptr ensures lifetime management.
 *
 * Contains:
 *   - all sequences,
 *   - all adduct/ion names,
 *   - all PTM names referenced by ModifiedPeptide.
 */
struct AnnotationResources {
    std::shared_ptr<std::vector<std::string>> sequences;
    std::shared_ptr<std::vector<std::string>> ions;
    std::shared_ptr<std::vector<std::string>> ptms;

    AnnotationResources()
        : sequences(std::make_shared<std::vector<std::string>>()),
          ions(std::make_shared<std::vector<std::string>>()),
          ptms(std::make_shared<std::vector<std::string>>()) {}
};

/**
* @brief Construct AnnotationInfo objects from PeptideEntry structures.
 *
 * Rebuilds:
 *   - sequence index
 *   - ModifiedPeptide (residue + terminal PTMs)
 *   - composition
 *   - decoy flag
 *
 */
struct AnnotationInfo {
    std::string* sequence;                     // points into resources->sequences
    std::shared_ptr<ModifiedPeptide> peptide;      // mods info
    std::string* adduct_name;                  // points into resources->ions
    charge_t charge;
    Composition total_composition;
    bool is_decoy;
    std::shared_ptr<AnnotationResources> resources; // keeps all strings alive

    AnnotationInfo(std::string* seq,
                   std::shared_ptr<ModifiedPeptide> pep,
                   std::string* adduct,
                   charge_t c,
                   const Composition& comp,
                   bool decoy,
                   std::shared_ptr<AnnotationResources> res)
        : sequence(seq),
          peptide(std::move(pep)),
          adduct_name(adduct),
          charge(c),
          total_composition(comp),
          is_decoy(decoy),
          resources(std::move(res)) {}
};

/**
 * @brief Writes annotation information to a CSV file.
 *
 * Columns include:
 *   - sequence
 *   - modifications
 *   - adduct
 *   - charge
 *   - monoisotopic mass
 *   - expected isotopic value
 *   - decoy flag
 */
void write_annotations_to_csv(
    const std::string& filename,
    const std::vector<AnnotationInfo>& annotations);


/**
 * @brief Generates AnnotationInfo objects for all peptide entries.
 *
 * Uses shared string resources to avoid duplication. Resolves:
 *   - peptide sequences,
 *   - adduct names,
 *   - modifications from ModifiedPeptide.
 *
 * @param fasta_records Parent sequences.
 * @param entries Peptide entries to annotate.
 * @return Vector of AnnotationInfo.
 */
std::vector<AnnotationInfo> generate_annotations(
    const std::vector<FastaRecord>& fasta_records,
    const std::vector<PeptideEntry>& entries);

/**
 * @brief Produces a compact annotation representation.
 *
 * Returns indices mapping peptide entries to FASTA references without
 * constructing full AnnotationInfo objects.
 *
 * @param fasta_records FASTA records.
 * @param entries Peptide entries.
 * @return Vector of annotation indices.
 */
std::vector<size_t> generate_better_annotations(
    const std::vector<FastaRecord>& fasta_records,
    const std::vector<PeptideEntry>& entries);

} // namespace wasabi::iso
