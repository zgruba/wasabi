#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <cctype>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>

#include <omp.h>
#include <thread>

#include "iso.hpp"

#include "amino_acids.hpp"
#include "wasabi_db.hpp"

namespace wasabi::ptm {
struct PTM;
}

template<typename T>
void write_int(std::ofstream& out, const T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template<typename T>
__always_inline T read_int(std::istream& in) {
    T val;
    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
}

template <typename T>
__always_inline void write_vector(std::ofstream& out, const std::vector<T>& v)
{
    const size_t n = v.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    if (n > 0)
        out.write(reinterpret_cast<const char*>(v.data()), n * sizeof(T));
}

template <typename T>
__always_inline std::vector<T> read_vector(std::ifstream& in)
{
    size_t n = 0;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    std::vector<T> v(n);
    if (n > 0)
        in.read(reinterpret_cast<char*>(v.data()), n * sizeof(T));
    return v;
}

template <typename T>
__always_inline void write_vector_raw(std::ofstream& out, const std::vector<T>& v)
{
    if (!v.empty())
        out.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
}

template <typename T>
__always_inline void read_vector_raw(std::ifstream& in, std::vector<T>& v)
{
    if (!v.empty())
        in.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(T));
}

namespace wasabi::iso {
struct ModifiedPeptide;

namespace fs = std::filesystem;

// ---------------- Ionisation modes ----------------
const std::map<std::string, wasabi::iso::IonisationMode> kIonisationModes = {
    // {"electron_loss_or_gain", {+1, 0.0}},
    {"protonation", {+1, 1.007276}},
    // {"deprotonation", {-1, -1.007276}},
    // {"Na_adduct", {+1, 22.989218}},
    // {"K_adduct", {+1, 38.963158}},
    // {"NH4_adduct", {+1, 18.033823}},
    // {"Li_adduct", {+1, 7.016004}},
    // {"Ca_adduct", {+2, 39.962591}},
    // {"Mg_adduct", {+2, 23.985042}},
    // {"Cl_adduct", {-1, 34.968853}},
    // {"formate_adduct", {-1, 44.998201}},
    // {"acetate_adduct", {-1, 59.013851}}
};

// -----------------------------------------
// Compute composition of a peptide sequence
// -----------------------------------------
inline Composition peptide_to_composition(std::string_view seq) {
    Composition total;
    for (const char aa : seq) {
        if (aa < 'A' || aa > 'Z')
            throw std::runtime_error(std::string("Invalid amino acid: ") + aa);

        total += aa::GetAminoAcidComposition(aa);
    }
    // Add water
    total += Composition(0,2,0,1,0,0,0);
    return total;
}

// -----------------------------------------
// Compute composition of a modified peptide
// -----------------------------------------
Composition compute_total_composition(const ModifiedPeptide& peptide)
{
    Composition total;

    // --- Add sequence composition ---
    for (const char aa : *peptide.sequence)
        total += aa::GetAminoAcidComposition(aa);  // existing function mapping AA to composition

    // add water
    total += Composition(0,2,0,1,0,0,0);

    // --- Add residue modifications (count-based) ---
    for (uint8_t i = 0; i < peptide.residue_mods_size; i++) {
        const auto& [count, ptmIdx] = peptide.residue_mods_ptr[i];
        total += ptm::common_ptms_k[ptmIdx].delta * count;
    }

    // --- Add terminal modifications ---
    if (peptide.has_n_term()) {
       total += peptide.get_n_term().delta;
    }

    if (peptide.has_c_term()) {
        total += peptide.get_c_term().delta;
    }

    return total;
}

// ----------------------------------------
// Get formula string (like "C12H22N4O6S1")
// ----------------------------------------
inline std::string compute_total_formula(const ModifiedPeptide& pep) {
    return compute_total_composition(pep).to_string();
}

// ----------------------------------------
// Generate modified peptides
// ----------------------------------------
std::vector<ModifiedPeptide> generate_modified_peptides_site(
    const std::string& seq,
    const uint8_t max_variable_mods)
{
    using namespace wasabi::ptm;
    std::vector<ModifiedPeptide> results;

    // --- Fixed PTMs ---
    std::vector<std::pair<uint8_t, common_ptms_size_t>> fixed_counts;
    for (const auto& ptm : common_ptms_k) {
        if (ptm.type == PtmType::Fixed && ptm.targetType == PtmTarget::Residue) {
            size_t count = std::count(seq.begin(), seq.end(), ptm.targetResidue);
            if (count > 0) {
                fixed_counts.emplace_back(count, ptm.common_index);
            }
        }
    }

    // --- Variable PTMs by type (residue) ---
    std::vector<std::pair<uint8_t, common_ptms_size_t>> variable_ptms; // (max_count, PTM)
    for (const auto& ptm : common_ptms_k) {
        if (ptm.type == PtmType::Variable) {
            if (ptm.targetType == PtmTarget::Residue) {
                size_t count = std::count(seq.begin(), seq.end(), ptm.targetResidue);
                if (count > 0) {
                    variable_ptms.emplace_back(count, ptm.common_index);
                }
            } else if (ptm.targetType == PtmTarget::ProteinNTerm || ptm.targetType == PtmTarget::ProteinCTerm) {
                // Treat terminal PTMs as count 1, since we can either apply or not
                variable_ptms.emplace_back(1, ptm.common_index);
            }
        }
    }

    // --- Recursive generation ---
    std::vector<std::pair<uint8_t, common_ptms_size_t>> current;

    std::function<void(uint8_t, uint8_t)> recurse;
    recurse = [&](const uint8_t idx, const uint8_t used) {
        if (used > max_variable_mods) return;

        if (idx == variable_ptms.size()) {
            // Merge fixed + variable counts
            std::vector<std::pair<uint8_t, common_ptms_size_t>> mods = fixed_counts;
            common_ptms_size_t c_term_mod_idx = no_op_ptm_idx;
            common_ptms_size_t n_term_mod_idx = no_op_ptm_idx;

            for (const auto& [_, ptmIdx] : current) {
                switch (common_ptms_k[ptmIdx].targetType) {
                    case PtmTarget::Residue: {
                        bool found = false;
                        for (auto& [count, p] : mods) {
                            if (p == ptmIdx) { count++; found = true; break; }
                        }
                        if (!found) mods.emplace_back(1, ptmIdx);
                        break;
                    } case PtmTarget::ProteinNTerm: {
                        n_term_mod_idx = ptmIdx;
                        break;
                    } case PtmTarget::ProteinCTerm: {
                        c_term_mod_idx = ptmIdx;
                        break;
                    }
                }
            }

            results.emplace_back(&seq, std::move(mods), c_term_mod_idx, n_term_mod_idx);
            return;
        }

        const auto max_count = variable_ptms[idx].first;
        const auto ptmIdx = variable_ptms[idx].second;

        // Try all counts from 0 to min(max_count, remaining slots)
        const auto remaining = static_cast<uint8_t>(max_variable_mods - used);
        const auto limit = std::min(max_count, remaining);
        for (uint8_t c = 0; c <= limit; ++c) {
            for (uint8_t i = 0; i < c; ++i) {
                current.emplace_back(0, ptmIdx);
            }
            recurse(idx + 1, used + c);
            for (uint8_t i = 0; i < c; ++i) {
                current.pop_back();
            }
        }
    };

    recurse(0, 0);

    return results;
}

// ----------------------------------------
// Generate entries from single record
// ----------------------------------------
void generate_spectra_from_single_record(
    const FastaRecord& record,
    const std::vector<multi_t>& multipliers,
    const size_t max_variable_mods,
    std::vector<PeptideEntry>& spectra)
{
    // Generate all modified peptides for this record
    auto modified_peptides = generate_modified_peptides_site(record.sequence, max_variable_mods);

    for (const auto& modpep : modified_peptides) {
        const Composition total = compute_total_composition(modpep);
        auto mod_pep_ptr = boost::intrusive_ptr(new ModifiedPeptide(modpep));

        for (const multi_t mult : multipliers) {
            for (const auto& [adduct_name, ion_info] : kIonisationModes) {
                // Copy-based shared_ptr for PeptideEntry
                spectra.emplace_back(
                    &record.id,
                    mod_pep_ptr,
                    static_cast<charge_t>(mult * ion_info.charge),
                    &adduct_name,
                    total,
                    record.is_decoy
                );
            }
        }
    }

    // Free memory used by the temporary modified_peptides vector
    std::vector<ModifiedPeptide>().swap(modified_peptides);

}

// ----------------------------------------
// Generate entries from fasta file
// ----------------------------------------
std::vector<PeptideEntry> generate_entries_from_fasta_records(
    const std::vector<FastaRecord>& records,
    const std::vector<multi_t>& multipliers,
    size_t max_variable_mods,
    int threads)
{
    if (records.empty()) return {};

    const int num_threads = std::max(1, threads);
    std::vector<std::vector<PeptideEntry>> local_spectra(static_cast<size_t>(num_threads));

#pragma omp parallel num_threads(num_threads) default(none) shared(local_spectra, records, multipliers, max_variable_mods)
    {
        const int tid = omp_get_thread_num();
        auto& local = local_spectra[static_cast<size_t>(tid)];


#pragma omp for schedule(static)
        for (const auto & record : records) {
            generate_spectra_from_single_record(record, multipliers, max_variable_mods, local);
        }
    }

    std::vector<PeptideEntry> spectra;

    // Merge all thread-local vectors into final vector and free memory
    for (auto& v : local_spectra) {
        spectra.insert(spectra.end(), v.begin(), v.end());
        std::vector<PeptideEntry>().swap(v);
    }
    spectra.shrink_to_fit();

    return spectra;
}

// ----------------------------------------
// Compute mz from mass and total charge
// ----------------------------------------
double compute_mz(double mass, const IonisationMode& ion_info, charge_t total_charge)
{
    // total_charge is already multiplied (e.g. +2, +3)
    // total_charge % ion_info.charge == 0
    const double total_shift = ion_info.mass_shift * (total_charge / ion_info.charge);

    if (total_charge != 0)
        return (mass + total_shift) / std::abs(total_charge);
    else
        return mass;
}

// ----------------------------------------
// Writes peptide entries to a binary file
// ----------------------------------------
void write_binary_info(
    const std::string& filename,
    const std::vector<FastaRecord>& fasta_records,
    const std::vector<PeptideEntry>& entries)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file for writing");

    // --- Precompute PTM indices ---
    std::vector<ptm_id_t> ptm_index_map(ptm::common_ptms_size_k);
    for (ptm_id_t i = 0; i < static_cast<ptm_id_t>(ptm::common_ptms_size_k); ++i)
        ptm_index_map[i] = i;

    // --- Precompute sequence indices ---
    std::unordered_map<std::string_view, seq_id_t> seq_index_map;
    for (seq_id_t i = 0; i < static_cast<seq_id_t>(fasta_records.size()); ++i)
        seq_index_map[fasta_records[i].sequence] = i;

    // --- Sequences ---
    for (const auto& rec : fasta_records)
        out.write(rec.sequence.data(), rec.sequence.size()), out.put('\0');
    out.put('\0');

    // --- PTMs ---
    for (const auto& ptm : ptm::common_ptms_k)
        out.write(ptm.name, strlen(ptm.name)), out.put('\0');
    out.put('\0');

    // --- Ions ---
    for (const auto& [name, _] : kIonisationModes)
        out.write(name.c_str(), name.size()), out.put('\0');
    out.put('\0');

    // --- Peptides ---
    for (const auto& entry : entries) {
        seq_id_t seq_idx = seq_index_map.at(*entry.peptide->sequence);
        write_int<seq_id_t>(out, seq_idx);
        write_int<ion_id_t>(out, 0); // placeholder
        write_int<charge_t>(out, static_cast<charge_t>(entry.charge));

        // --- Aggregate residue PTM counts in a vector ---
        std::vector<uint16_t> ptm_counts(ptm::common_ptms_size_k, 0);
        for (uint8_t i = 0; i < entry.peptide->residue_mods_size; i++) {
            const auto& [count, ptmIdx] = entry.peptide->residue_mods_ptr[i];
            ptm_counts[ptmIdx] += count;
        }
        // Count non-zero PTMs
        uint16_t non_zero_count = 0;
        for (auto c : ptm_counts) if (c) ++non_zero_count;
        write_int<uint16_t>(out, non_zero_count);

        for (ptm_id_t i = 0; i < static_cast<ptm_id_t>(ptm_counts.size()); ++i) {
            if (ptm_counts[i] == 0) continue;
            write_int<uint16_t>(out, ptm_counts[i]);
            write_int<ptm_id_t>(out, i);
        }

        // --- Terminal mods ---
        write_int<uint16_t>(out, entry.peptide->terminal_mods_count());
        if (entry.peptide->has_n_term()) {
            write_int<uint16_t>(out, entry.peptide->n_term_mod_idx);
        }
        if (entry.peptide->has_c_term()) {
            write_int<uint16_t>(out, entry.peptide->c_term_mod_idx);
        }

        // --- Decoy flag ---
        out.put(static_cast<char>(entry.is_decoy ? 1 : 0));
    }
}

// ----------------------------------------------
// Read annotations written by write_binary_info
// ----------------------------------------------
std::vector<AnnotationInfo> read_from_binary_info(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file for reading");

    auto resources = std::make_shared<AnnotationResources>();
    std::string str;

    // --- Helper to read strings into shared_ptr<vector<string>> ---
    auto read_strings = [&](std::shared_ptr<std::vector<std::string>>& container) {
        while (std::getline(in, str, '\0') && !str.empty())
            container->push_back(str);
    };

    read_strings(resources->sequences);
    read_strings(resources->ptms);
    read_strings(resources->ions);

    std::vector<AnnotationInfo> annots;

    while (in && !in.eof()) {
        seq_id_t seq_idx = read_int<seq_id_t>(in);
        ion_id_t ion_idx = read_int<ion_id_t>(in);
        charge_t charge = read_int<charge_t>(in);

        // --- Residue mods ---
        uint16_t mod_count = read_int<uint16_t>(in);
        std::vector<std::pair<uint8_t, ptm::common_ptms_size_t>> residue_mods;
        residue_mods.reserve(mod_count);
        for (uint16_t i = 0; i < mod_count; ++i) {
            size_t count = read_int<uint16_t>(in);
            ptm_id_t idx = read_int<ptm_id_t>(in);
            residue_mods.emplace_back(count, idx);
        }

        // --- Terminal mods ---
        ptm::common_ptms_size_t c_term_mod_idx = ptm::no_op_ptm_idx;
        ptm::common_ptms_size_t n_term_mod_idx = ptm::no_op_ptm_idx;
        uint16_t term_count = read_int<uint16_t>(in);
        for (uint16_t i = 0; i < term_count; ++i) {
            switch (ptm_id_t idx = read_int<ptm_id_t>(in); ptm::common_ptms_k[idx].targetType) {
                case ptm::PtmTarget::ProteinCTerm: {
                    c_term_mod_idx = idx;
                    break;
                } case ptm::PtmTarget::ProteinNTerm: {
                    n_term_mod_idx = idx;
                    break;
                } default: {
                    throw std::runtime_error("Incorrect ptm target type");
                }
            }
        }

        // --- Decoy flag ---
        char flag = 0;
        if (!in.get(flag)) break;
        bool is_decoy = (flag != 0);

        auto pep = std::make_shared<ModifiedPeptide>(
            ModifiedPeptide(
                &resources->sequences->at(seq_idx),
                std::move(residue_mods),
                c_term_mod_idx,
                n_term_mod_idx
            )
        );

        Composition comp = compute_total_composition(*pep);

        annots.emplace_back(
            &resources->sequences->at(seq_idx),
            pep,
            &resources->ions->at(ion_idx),
            charge,
            comp,
            is_decoy,
            resources
        );

    }

    return annots;
}


// ------------------------------------
// Generate only the isotopic info
// ------------------------------------
IsotopicInfo generate_isotopic_info(const PeptideEntry& pep, double threshold, size_t id)
{
    using namespace IsoSpec;

    // Generate formula string from Composition
    std::string formula = pep.total_composition.to_string();

    Iso iso(formula.c_str());
    IsoThresholdGenerator gen(std::move(iso), (1 - threshold) / 1e2);

    std::vector<double> masses;
    std::vector<double> probs;

    double total_prob = 0.0;
    while (gen.advanceToNextConfiguration() && total_prob < threshold)
    {
        double p = gen.prob();
        masses.push_back(gen.mass());
        probs.push_back(p);
        total_prob += p;
    }

    if (total_prob == 0.0)
    {
        return {id, 0.0, 0, pep.adduct_name, pep.charge};
    }

    for (auto& p : probs)
        p /= total_prob;

    // Weighted average mass
    double expected_mass = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        expected_mass += masses[i] * probs[i];

    // Compute m/z based on ionization mode and charge
    const auto& ion_mode = kIonisationModes.at(*pep.adduct_name);
    double expected_mz = compute_mz(expected_mass, ion_mode, pep.charge);

    return {id, expected_mz, static_cast<dist_size_t>(masses.size()), pep.adduct_name, pep.charge};
}

// ----------------------------------------------------
// Write full isotopic distributions to a binary file.
// ----------------------------------------------------
void write_binary_distribution(const std::string& filename,
                               const std::vector<PeptideEntry>& entries,
                               const std::vector<IsotopicInfo>& isotopic_infos,
                               double threshold)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Cannot open file for writing: " + filename);

    const auto n_entries = static_cast<entry_id_t>(isotopic_infos.size());
    write_int(out, n_entries);

    std::vector<HeaderRecord> headers;
    headers.reserve(n_entries);

    const size_t header_section_size = n_entries * sizeof(HeaderRecord);
    auto current_offset = sizeof(n_entries) + header_section_size;

    for (const auto& info : isotopic_infos) {
        if (info.distribution_size == 0)
            throw std::runtime_error("Cannot write isotopic info with zero distribution size");

        headers.emplace_back(info.id, info.expected, info.distribution_size, current_offset);
        current_offset += info.distribution_size * (sizeof(mz_t) + sizeof(intensity_t));
    }

    // --- Write headers ---
    write_vector_raw(out, headers);

    // --- Write data ---
    std::vector<mz_t> mzs;
    std::vector<intensity_t> probs;

    for (const auto& info : isotopic_infos)
    {
        if (info.distribution_size == 0)
            continue;

        const auto& e = entries[info.id];
        std::string formula = e.total_composition.to_string(); // updated line

        IsoSpec::Iso iso(formula.c_str());
        IsoSpec::IsoThresholdGenerator gen(std::move(iso), (1 - threshold) / 1e2);

        mzs.resize(info.distribution_size);
        probs.resize(info.distribution_size);

        double total_prob = 0.0;
        size_t pos = 0;

        while (gen.advanceToNextConfiguration() && total_prob < threshold)
        {
            double prob = gen.prob();
            double mass = gen.mass();
            double mz = compute_mz(mass, kIonisationModes.at(*info.adduct_name), info.total_charge);

            mzs[pos] = mz;
            probs[pos] = prob;
            total_prob += prob;
            ++pos;
        }

        // Trim unused space
        mzs.resize(pos);
        probs.resize(pos);

        // Sort isotopes by increasing m/z
        std::vector<size_t> indices(pos);
        std::iota(indices.begin(), indices.end(), 0);

        std::sort(indices.begin(), indices.end(),
                  [&](size_t a, size_t b) { return mzs[a] < mzs[b]; });

        std::vector<mz_t> sorted_mzs(pos);
        std::vector<intensity_t> sorted_probs(pos);
        for (size_t i = 0; i < pos; ++i) {
            sorted_mzs[i] = mzs[indices[i]];
            sorted_probs[i] = probs[indices[i]];
        }

        if (total_prob > 0.0)
            for (auto& p : sorted_probs)
                p /= total_prob;

        // Write both sorted vectors
        write_vector_raw(out, sorted_mzs);
        write_vector_raw(out, sorted_probs);
    }
}

// ----------------------------------------------------
// Write full isotopic distributions to a binary file.
// ----------------------------------------------------
void write_annotations_to_csv(
    const std::string& filename,
    const std::vector<AnnotationInfo>& annotations)
{
    std::ofstream csv(filename);
    if (!csv) throw std::runtime_error("Cannot open CSV for writing");

    csv << "SEQUENCE,ADDUCT,CHARGE,MODIFICATIONS,TOTAL_COMPOSITION,IS_DECOY\n";

    for (const auto& a : annotations) {
        csv << *a.sequence << ",";
        csv << *a.adduct_name << ",";
        csv << std::to_string(a.charge) << ",\"";

        if (a.peptide) {
            bool first = true;
            for (uint8_t i = 0; i < a.peptide->residue_mods_size; i++) {
                const auto& [pos, ptmIdx] = a.peptide->residue_mods_ptr[i];
                if (!first) csv << ";";
                csv << ptm::common_ptms_k[ptmIdx].name << "@" << std::to_string(pos);
                first = false;
            }
            if (a.peptide->has_c_term()) {
                if (!first) csv << ";";
                csv << a.peptide->get_c_term().name << "@C-term";
                first = false;
            }
            if (a.peptide->has_n_term()) {
                if (!first) csv << ";";
                csv << a.peptide->get_n_term().name << "@N-term";
            }
        }

        csv << "\",";
        csv << a.total_composition.to_string() << ",";
        csv << (a.is_decoy ? "true" : "false") << "\n";
    }
}

// ------------------------------------------------
// Sorts peptides by their expected isotopic mass
// ------------------------------------------------
std::vector<IsotopicInfo>
sort_peptides_by_expected_value(const std::vector<PeptideEntry>& peptides, const double threshold) {
    std::vector<IsotopicInfo> infos;
    infos.reserve(peptides.size());

    for (size_t i = 0; i < peptides.size(); ++i) {
        infos.push_back(generate_isotopic_info(peptides[i], threshold, i));
    }

    std::sort(infos.begin(), infos.end(),
              [](const IsotopicInfo& a, const IsotopicInfo& b) {
                  return a.expected < b.expected;
              });

    return infos;
}

// -------------------------------------------
// Read binary file with entries distribution
// -------------------------------------------
std::vector<std::tuple<size_t, mz_t, std::vector<mz_t>, std::vector<intensity_t>>>
read_binary_distribution(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }

    // --- Number of entries ---
    auto n_entries = read_int<entry_id_t>(in);

    // --- Read header table ---
    std::vector<HeaderRecord> headers(n_entries);
    read_vector_raw(in, headers);

    // --- Read each data block ---
    std::vector<std::tuple<size_t, mz_t, std::vector<mz_t>, std::vector<intensity_t>>> results;
    results.reserve(n_entries);

    for (const auto& h : headers) {
        in.seekg(h.offset, std::ios::beg);

        std::vector<mz_t> mzs(h.size);
        std::vector<intensity_t> probs(h.size);

        read_vector_raw(in, mzs);
        read_vector_raw(in, probs);

        results.emplace_back(h.id, h.expected, std::move(mzs), std::move(probs));
    }

    return results;
}

template <typename T>
void apply_permutation_in_place(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i])
        {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

// -----------------------------------------------
// Reads full isotopic distributions from binary
// -----------------------------------------------
void generate_isotopic_distribution(const PeptideEntry &pep, const double threshold,
                                    std::vector<double> &mzs, std::vector<double> &probs) {
    IsoSpec::Iso iso(pep.total_composition.to_string());
    IsoSpec::IsoThresholdGenerator gen(std::move(iso), (1 - threshold) / 1e2);
    const auto &ion_mode = kIonisationModes.at(*pep.adduct_name);
    const auto charge = pep.charge;

    double total_prob = 0.0;
    while (gen.advanceToNextConfiguration() && total_prob < threshold) {
        double p = gen.prob();
        mzs.push_back(compute_mz(gen.mass(), ion_mode, charge));
        probs.push_back(p);
        total_prob += p;
    }

    for (auto &p: probs) {
        p /= total_prob;
    }

    auto indexes = std::vector<size_t>(mzs.size(), 0);
    std::iota(indexes.begin(), indexes.end(), 0);
    std::sort(indexes.begin(), indexes.end(), [&](const size_t a, const size_t b) { return mzs[a] < mzs[b]; });
    apply_permutation_in_place( mzs, indexes);
    apply_permutation_in_place( probs, indexes);
}

// ------------------------------------------------
// Generates AnnotationInfo for all peptide entries
// ------------------------------------------------
std::vector<AnnotationInfo> generate_annotations(const std::vector<FastaRecord>& fasta_records,
                     const std::vector<PeptideEntry>& entries)
{
    auto resources = std::make_shared<AnnotationResources>();

    for (const auto& rec : fasta_records)
        resources->sequences->push_back(rec.sequence);

    for (const auto& p : ptm::common_ptms_k)
        resources->ptms->emplace_back(p.name);

    for (const auto& [name, _] : kIonisationModes)
        resources->ions->push_back(name);

    std::vector<AnnotationInfo> annots;
    annots.reserve(entries.size());

    std::unordered_map<std::string_view, seq_id_t> seq_index;
    seq_index.reserve(fasta_records.size());

    for (seq_id_t i = 0; i < fasta_records.size(); ++i)
        seq_index[fasta_records[i].sequence] = i;

    for (const auto& entry : entries)
    {
        const auto& seq_sv = *entry.peptide->sequence;
        auto it = seq_index.find(seq_sv);
        if (it == seq_index.end())
            throw std::runtime_error("PeptideEntry sequence not present in FASTA table");

        const seq_id_t seq_idx = it->second;

        std::vector residue_mods(
            entry.peptide->residue_mods_ptr,
            entry.peptide->residue_mods_ptr + entry.peptide->residue_mods_size
        );

        auto pep = std::make_shared<ModifiedPeptide>(
            &resources->sequences->at(seq_idx),
            std::move(residue_mods),
            entry.peptide->c_term_mod_idx,
            entry.peptide->n_term_mod_idx
        );

        Composition comp = compute_total_composition(*pep);

        ion_id_t ion_idx = 0;

        annots.emplace_back(
            &resources->sequences->at(seq_idx),  // pointer to sequence
            pep,                                  // modified peptide
            &resources->ions->at(ion_idx),        // ion name
            static_cast<charge_t>(entry.charge),  // charge
            comp,                                 // chemical composition
            entry.is_decoy,                       // decoy flag
            resources                             // shared annotation resources
        );
    }

    return annots;
}

// -------------------------------------------
// Produces a compact annotation representation
// -------------------------------------------
std::vector<size_t> generate_better_annotations(const std::vector<FastaRecord>& fasta_records,
                     const std::vector<PeptideEntry>& entries) {
    std::unordered_map<std::string_view, seq_id_t> seq_index;
    for (seq_id_t i = 0; i < fasta_records.size(); ++i) {
        seq_index[fasta_records[i].sequence] = i;
    }
    std::vector<size_t> annots;
    annots.reserve(entries.size());
    for (const auto& entry : entries) {
        annots.push_back(seq_index[*entry.peptide->sequence]);
    }
    return annots;
}

}
