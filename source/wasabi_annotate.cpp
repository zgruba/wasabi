#include <ranges>
#include <string>
#include <vector>
#include <iomanip>
#include <limits>

#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>

#include "iso.hpp"
#include "parsing.hpp"

namespace po = boost::program_options;

using PeptideEntry = wasabi::iso::PeptideEntry;
using PeptideEntries = std::vector<PeptideEntry>;
using compressed_t = std::vector<std::tuple<const PeptideEntry*, size_t, size_t>>;
using db_index_t = std::vector<std::pair<double, std::size_t>>;
using database_t = wasabi::db::StreamQWassersteinDatabase<double, double>;
using queries_t = std::vector<wasabi::mass::CompactFeatureSpec>;
using search_result_t = std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<double>>>; // Indices, Distances.

constexpr bool bin_mz_k = false;
constexpr double bin_width_k = 0.01;
constexpr int missed_cleavages_k = 2;
constexpr int sequence_min_len_k = 7;
constexpr int sequence_max_len_k = 51;
constexpr bool allow_semitryptic_k = false;
constexpr int max_variable_modifications_k = 3;
constexpr double isotopic_threshold_k = 0.999;

struct Modification {
    std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>* residue_mods_ptr = nullptr;
    wasabi::ptm::common_ptms_size_t residue_mods_size = 0;
    wasabi::ptm::common_ptms_size_t c_term_mod_idx = wasabi::ptm::no_op_ptm_idx;
    wasabi::ptm::common_ptms_size_t n_term_mod_idx =  wasabi::ptm::no_op_ptm_idx;

    explicit Modification(const std::string& mods_str);
    explicit Modification(const wasabi::iso::ModifiedPeptide& modified_peptide);
    Modification(const Modification& other) noexcept;
    Modification(const Modification&& other) noexcept;
    ~Modification();
};

struct Config {
    std::string feature_path; // necessary
    std::string mzml_path; // necessary
    std::string fasta_path; // necessary
    std::string rt_mapping_path; // optional; if use-rt-filtering set to true then necessary
    std::string output_path; // necessary
    std::vector<int8_t> multipliers; // necessary
    double queue_threshold = 0.15; // optional
    double search_threshold = 0.15; // optional
    double rt_threshold = 1750; // optional
    bool reduced_decoys = false; // optional
    bool with_decoys = false; // optional
    bool use_rt_filtering = false; // optional
    int threads = 1; // optional
};

void log_config(const Config& config);

void print_usage_and_exit(const po::options_description& all, int exit_code);

Config parse_args(int argc, char* argv[]);

std::pair<const wasabi::FastaRecords, PeptideEntries> generate_entries(
        const std::string& fasta_file, const std::vector<int8_t>& multipliers, bool with_decoys,
        bool reduced_decoys, int threads);

PeptideEntries reduce_decoys(PeptideEntries& entries);

void enrich_entries_with_retention_times(const std::string& file_with_rt, PeptideEntries& entries, int threads);

compressed_t compress_entries(const PeptideEntries& entries);

db_index_t create_database_index(const compressed_t& entries, int threads);

database_t create_database(const db_index_t& db_index, const compressed_t& entries,
                           double queue_threshold, double search_threshold);

queries_t load_queries(const std::string& feature_file, const std::string& mzml_file, int threads);

search_result_t answer_queries(const database_t& db,
                               const PeptideEntries& entries,
                               const queries_t& queries,
                               const compressed_t& compressed_entries,
                               double rt_threshold,
                               int threads);

void save_annotations(const wasabi::FastaRecords& records,
                      const PeptideEntries& entries,
                      const std::string& output_file,
                      const queries_t& queries,
                      const search_result_t& search_result);

int main(const int argc, char* argv[]) {
    const auto cfg = parse_args(argc, argv);
    log_config(cfg);

    auto [records, entries] = generate_entries(
        cfg.fasta_path, cfg.multipliers, cfg.with_decoys, cfg.reduced_decoys, cfg.threads);

    if (cfg.use_rt_filtering) {
        enrich_entries_with_retention_times(cfg.rt_mapping_path, entries, cfg.threads);
    }

    const auto compressed_entries = compress_entries(entries);

    const auto db_index = create_database_index(compressed_entries, cfg.threads);

    const auto database = create_database(db_index, compressed_entries, cfg.queue_threshold, cfg.search_threshold);

    const auto queries = load_queries(cfg.feature_path, cfg.mzml_path, cfg.threads);

    const auto search_result = answer_queries(
        database, entries, queries, compressed_entries, cfg.rt_threshold, cfg.threads);

    save_annotations(records, entries, cfg.output_path, queries, search_result);

    BOOST_LOG_TRIVIAL(info) << "Finished annotating.";
}

boost::log::basic_record_ostream<char>& operator<<(boost::log::basic_record_ostream<char>& os, const std::vector<int8_t>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << static_cast<int>(v[i]);
        if (i + 1 < v.size()) os << ", ";
    }
    os << "]";
    return os;
}

void log_config(const Config& config) {
    BOOST_LOG_TRIVIAL(info) << "Running wasabi annotator with parameters:"
    << " threads=" << config.threads
    << ", with_decoys=" << config.with_decoys
    << ", reduced_decoys=" << config.reduced_decoys
    << ", multipliers=" << config.multipliers
    << ", queue_threshold=" << config.queue_threshold
    << ", search_threshold=" << config.search_threshold
    << ", rt_threshold=" << config.rt_threshold
    << ", use_rt_filtering=" << config.use_rt_filtering
    << ".";
    BOOST_LOG_TRIVIAL(info) << "Input files are:";
    BOOST_LOG_TRIVIAL(info) << "  feature_file=" << config.feature_path;
    BOOST_LOG_TRIVIAL(info) << "  mzml_file=" << config.mzml_path;
    BOOST_LOG_TRIVIAL(info) << "  fasta_file=" << config.fasta_path;
    if (config.use_rt_filtering) {
        BOOST_LOG_TRIVIAL(info) << "  rt_mapping_file=" << config.rt_mapping_path;
    }
    BOOST_LOG_TRIVIAL(info) << "Output file is: " << config.output_path;
}

void print_usage_and_exit(const po::options_description& all, int exit_code) {
    std::cout << "Usage:\n"
              << "  wasabi_annotate <features> <mzML> <fasta> <output> "
              << "--multipliers <list> [options]\n\n"
              << all << "\n";
    std::exit(exit_code);
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;

    // Named options
    po::options_description named_opts("Options");
    named_opts.add_options()
        ("help,h", "Show help message")
        ("threads,t", po::value<int>()->default_value(1),
            "Number of threads (default: 1)")
        ("multipliers,m", po::value<std::vector<int>>()->multitoken(),
            "Multipliers (space-separated list, required)")
        ("queue-threshold", po::value<double>()->default_value(0.15),
            "Queue threshold (default: 0.15)")
        ("search-threshold", po::value<double>()->default_value(0.15),
            "Search threshold (default: 0.15)")
        ("rt-threshold", po::value<double>()->default_value(1750),
            "RT threshold (default: 1750)")
        ("reduced-decoys", po::bool_switch()->default_value(false),
            "Use reduced decoys")
        ("with-decoys,d", po::bool_switch()->default_value(false),
            "Include decoys")
        ("use-rt-filtering", po::bool_switch()->default_value(false),
            "Enable RT filtering")
        ("rt-mapping", po::value<std::string>(),
            "RT mapping file (required if --use-rt-filtering is set)");

    // Positional arguments
    po::options_description positional_opts("Positional arguments");
    positional_opts.add_options()
        ("feature_path", po::value<std::string>(), "Feature file (required)")
        ("mzml_path", po::value<std::string>(), "mzML file (required)")
        ("fasta_path", po::value<std::string>(), "FASTA file (required)")
        ("output_path", po::value<std::string>(), "Output path (required)");

    po::options_description all;
    all.add(named_opts).add(positional_opts);

    po::positional_options_description pos;
    pos.add("feature_path", 1);
    pos.add("mzml_path", 1);
    pos.add("fasta_path", 1);
    pos.add("output_path", 1);

    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all)
                      .positional(pos)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << "\n\n";
        print_usage_and_exit(all, 1);
    }

    if (vm.contains("help")) {
        print_usage_and_exit(all, 0);
    }

    // Validate required positional arguments
    const char* required_positional[] = {
        "feature_path",
        "mzml_path",
        "fasta_path",
        "output_path"
    };

    for (const char* key : required_positional) {
        if (!vm.contains(key)) {
            std::cerr << "Error: missing required positional argument: "
                      << key << "\n\n";
            print_usage_and_exit(all, 1);
        }
    }

    // Validate required named arguments
    if (!vm.contains("multipliers")) {
        std::cerr << "Error: --multipliers is required.\n\n";
        print_usage_and_exit(all, 1);
    }

    // Assign values
    cfg.feature_path = vm["feature_path"].as<std::string>();
    cfg.mzml_path = vm["mzml_path"].as<std::string>();
    cfg.fasta_path = vm["fasta_path"].as<std::string>();
    cfg.output_path = vm["output_path"].as<std::string>();

    cfg.queue_threshold = vm["queue-threshold"].as<double>();
    cfg.search_threshold = vm["search-threshold"].as<double>();
    cfg.rt_threshold = vm["rt-threshold"].as<double>();
    cfg.reduced_decoys = vm["reduced-decoys"].as<bool>();
    cfg.with_decoys = vm["with-decoys"].as<bool>();
    cfg.use_rt_filtering = vm["use-rt-filtering"].as<bool>();
    cfg.threads = vm["threads"].as<int>();

    // Multipliers conversion to int8_t
    const auto& mults = vm["multipliers"].as<std::vector<int>>();
    if (mults.empty()) {
        std::cerr << "Error: multipliers list must not be empty.\n\n";
        print_usage_and_exit(all, 1);
    }

    for (int v : mults) {
        if (v < std::numeric_limits<int8_t>::min() ||
            v > std::numeric_limits<int8_t>::max()) {
            std::cerr << "Error: multiplier out of int8_t range: "
                      << v << "\n\n";
            print_usage_and_exit(all, 1);
        }
        cfg.multipliers.push_back(static_cast<int8_t>(v));
    }

    // Conditional requirement
    if (cfg.use_rt_filtering) {
        if (!vm.contains("rt-mapping")) {
            std::cerr << "Error: --rt-mapping is required when "
                         "--use-rt-filtering is enabled.\n\n";
            print_usage_and_exit(all, 1);
        }
        cfg.rt_mapping_path = vm["rt-mapping"].as<std::string>();
    }

    // Basic validation
    if (cfg.threads <= 0) {
        std::cerr << "Error: threads must be >= 1.\n\n";
        print_usage_and_exit(all, 1);
    }

    return cfg;
}

PeptideEntries reduce_decoys(PeptideEntries& entries) {
    BOOST_LOG_TRIVIAL(info) << "Removing decoys with existing non-decoy equal entries.";
    std::ranges::sort(entries, [](const PeptideEntry& a, const PeptideEntry& b) {
        if (a.total_composition.elements < b.total_composition.elements) return true;
        if (a.total_composition.elements > b.total_composition.elements) return false;
        if (a.charge < b.charge) return true;
        if (a.charge > b.charge) return false;
        return  (a.is_decoy < b.is_decoy);
    });
    const auto eq = [](const PeptideEntry& a, const PeptideEntry& b) {
        return a.total_composition.elements == b.total_composition.elements && a.charge == b.charge;
    };
    // As this function was written after rejecting researching semi-tryptic peptides it can be further improved.
    PeptideEntries result;
    result.push_back(entries[0]);
    bool chunk_has_non_decoy = !entries[0].is_decoy;
    for (size_t i = 1; i < entries.size(); ++i) {
        if (eq(entries[i], entries[i - 1])) {
            if (entries[i].is_decoy) {
                if (!chunk_has_non_decoy) {
                    result.push_back(entries[i]);
                }
            } else {
                result.push_back(entries[i]);
                chunk_has_non_decoy = true;
            }
        } else {
            result.push_back(entries[i]);
            chunk_has_non_decoy = !entries[i].is_decoy;
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Removed " << entries.size() - result.size() << " entries.";
    return result;
}

std::pair<const wasabi::FastaRecords, PeptideEntries> generate_entries(
        const std::string& fasta_file, const std::vector<int8_t>& multipliers, const bool with_decoys,
        const bool reduced_decoys, const int threads) {
    BOOST_LOG_TRIVIAL(info) << "Generating entries.";

    auto records = wasabi::fasta::read_fasta_with_trypsin(
        fasta_file, missed_cleavages_k, sequence_min_len_k, sequence_max_len_k,allow_semitryptic_k, with_decoys);

    BOOST_LOG_TRIVIAL(info) << "Read " << records.size() << " records from the fasta file.";

    auto entries = wasabi::iso::generate_entries_from_fasta_records(
        records, multipliers, max_variable_modifications_k, threads);

    BOOST_LOG_TRIVIAL(info) << "Generated " << entries.size() << " entries.";

    if (with_decoys && reduced_decoys) {
       entries = reduce_decoys(entries);
    }

    std::sort(entries.begin(), entries.end());

    return {std::move(records), std::move(entries)};
}

std::vector<std::string> split_str(const std::string& s, const char delim) {
    std::vector<std::string> parts;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        parts.push_back(item);
    }
    return parts;
}

Modification::Modification(const std::string& mods_str) {
    std::vector<std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>> residue_mods;

    for (const auto & mod :  split_str(mods_str, ';')) {
        if (auto terms = split_str(mod, '@'); terms.size() == 2) {
            auto [_, idx] = wasabi::ptm::decode_modification(terms[0]);
            if (terms[1][0] == 'N') {
                n_term_mod_idx = idx;
            } else {
                c_term_mod_idx = idx;
            }
        } else {
            residue_mods.emplace_back(wasabi::ptm::decode_modification(mod));
        }
    }

    residue_mods_size = residue_mods.size();
    residue_mods_ptr = new std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>[residue_mods_size];
    for (int i = 0; i < residue_mods_size; ++i) {
        residue_mods_ptr[i] = residue_mods[i];
    }
}

Modification::Modification(const wasabi::iso::ModifiedPeptide& modified_peptide) : residue_mods_size(modified_peptide.residue_mods_size),
            c_term_mod_idx(modified_peptide.c_term_mod_idx), n_term_mod_idx(modified_peptide.n_term_mod_idx) {
        residue_mods_ptr = new std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>[residue_mods_size];
        for (int i = 0; i < residue_mods_size; ++i) {
            residue_mods_ptr[i] = modified_peptide.residue_mods_ptr[i];
        }
    }

Modification::Modification(const Modification& other) noexcept : residue_mods_size(other.residue_mods_size),
            c_term_mod_idx(other.c_term_mod_idx), n_term_mod_idx(other.n_term_mod_idx) {
        if (other.residue_mods_ptr) {
            residue_mods_ptr = new std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>[residue_mods_size];
            for (int i = 0; i < residue_mods_size; ++i) {
                residue_mods_ptr[i] = other.residue_mods_ptr[i];
            }
        } else {
            residue_mods_ptr = nullptr;
        }
    }

Modification::Modification(const Modification&& other) noexcept :
        residue_mods_size(other.residue_mods_size),
        c_term_mod_idx(other.c_term_mod_idx), n_term_mod_idx(other.n_term_mod_idx) {
        if (other.residue_mods_ptr) {
            residue_mods_ptr = new std::pair<uint8_t, wasabi::ptm::common_ptms_size_t>[residue_mods_size];
            for (int i = 0; i < residue_mods_size; ++i) {
                residue_mods_ptr[i] = other.residue_mods_ptr[i];
            }
        } else {
            residue_mods_ptr = nullptr;
        }
    }

Modification::~Modification() {
        delete[] residue_mods_ptr;
    }

bool operator==(const Modification& lhs, const Modification& rhs) {
    if (lhs.c_term_mod_idx == rhs.c_term_mod_idx && lhs.n_term_mod_idx == rhs.n_term_mod_idx &&
        lhs.residue_mods_size == rhs.residue_mods_size) {
        for (wasabi::ptm::common_ptms_size_t i = 0; i < lhs.residue_mods_size; ++i) {
            if (lhs.residue_mods_ptr[i] != rhs.residue_mods_ptr[i]) {
                return false;
            }
        }
        return true;
    }

    return false;
}

struct ModificationHash {
    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    size_t operator()(const Modification& m) const noexcept {
        uint64_t h = 0;

        h ^= mix(m.c_term_mod_idx + 0x100) ;
        h ^= mix((h << 1) ^ m.n_term_mod_idx + 0x200);
        h ^= mix((h << 1) ^ m.residue_mods_size);

        for (wasabi::ptm::common_ptms_size_t i = 0; i < m.residue_mods_size; ++i) {
            const auto&[fst, snd] = m.residue_mods_ptr[i];
            const uint64_t v = (static_cast<uint64_t>(fst) << 32) ^ static_cast<uint64_t>(snd);
            h ^= mix(v + i * 0x9e3779b9ULL);
        }

        return h;
    }
};

void enrich_entries_with_retention_times(const std::string& file_with_rt, PeptideEntries& entries, const int threads) {
    BOOST_LOG_TRIVIAL(info) << "Adding retention times to entries. Parsing the csv file.";

    std::unordered_map<std::string, std::unordered_map<Modification, double, ModificationHash>> rt_mapping;

    // Read csv into mapping
    {
        std::ifstream in_file_with_rt{file_with_rt};
        if (!in_file_with_rt) {
            throw std::runtime_error("Cannot open file for reading: " + file_with_rt);
        }
        std::string line;
        // Ignore first line:
        std::getline(in_file_with_rt, line);
        while (std::getline(in_file_with_rt, line)) {
            std::stringstream ss(line);
            std::string seq, mods, rt_str;

            std::getline(ss, seq, ',');
            std::getline(ss, mods, ',');
            std::getline(ss, rt_str, ',');
            double rt = std::stod(rt_str);

            rt_mapping[seq].emplace(mods, rt);
        }
    }

    size_t max_mods = 0;
    for (const auto &val: rt_mapping | std::views::values) {
        max_mods = std::max(max_mods, val.size());
    }
    BOOST_LOG_TRIVIAL(info) << "Read retention times of " << rt_mapping.size() << " unique sequences"
        << " with at most " << max_mods << " modifications each. Adding to entries.";

#pragma omp parallel num_threads(threads) default(none) shared(entries, rt_mapping)
    {
#pragma omp for schedule(static)
        for (auto& entry : entries) {
            if (auto seq_map = rt_mapping.find(*entry.peptide->sequence); seq_map != rt_mapping.end()) {
                if (auto rt = seq_map->second.find(Modification(*entry.peptide)); rt != seq_map->second.end()) {
                    entry.rt = rt->second;
                }
            }
        }
    }
}

compressed_t compress_entries(const PeptideEntries& entries) {
    BOOST_LOG_TRIVIAL(info) << "Compressing entries.";

    compressed_t compressed_entries;
    auto current_ptr = &entries[0];
    size_t count = 1, pos = 0;

    for (size_t i = 1; i < entries.size(); ++i) {
        if (entries[i] == *current_ptr) {
            ++count;
        } else {
            compressed_entries.emplace_back(current_ptr, pos, count);
            current_ptr = &entries[i];
            count = 1;
            pos = i;
        }
    }

    compressed_entries.emplace_back(current_ptr, pos, count);

    return compressed_entries;
}

db_index_t create_database_index(const compressed_t& entries, const int threads) {
    BOOST_LOG_TRIVIAL(info) << "Calculating expected values.";

    db_index_t db_index(entries.size());

#pragma omp parallel num_threads(threads) default(none) shared(db_index, entries)
    {
#pragma omp for schedule(static)
        for (size_t i = 0; i < entries.size(); ++i) {
            db_index[i] = {wasabi::iso::generate_isotopic_info(*std::get<0>(entries[i]), isotopic_threshold_k, i).expected, i};
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Expected values calculated. Sorting.";
    std::ranges::sort(db_index);

    return db_index;
}

database_t create_database(const db_index_t& db_index, const compressed_t& entries,
                           const double queue_threshold, const double search_threshold) {
    BOOST_LOG_TRIVIAL(info) << "Creating database.";

    auto dist_gen = [&](const size_t idx) -> wasabi::db::CompactDistribution<double,double> {
        thread_local std::vector<double> mzs_buf;
        thread_local std::vector<double> probs_buf;

        const auto& pep = *std::get<0>(entries[idx]);
        mzs_buf.clear();
        probs_buf.clear();

        wasabi::iso::generate_isotopic_distribution(pep, isotopic_threshold_k, mzs_buf, probs_buf);

        wasabi::db::CompactDistribution<double,double> dist{};
        dist.size = mzs_buf.size();
        dist.values = mzs_buf.data();
        dist.probabilities = probs_buf.data();

        return dist; // only valid while thread is using it
    };

    wasabi::db::StreamQWassersteinDatabase<double, double> db(db_index, dist_gen, queue_threshold, search_threshold);

    return db;
}

queries_t load_queries(const std::string& feature_file, const std::string& mzml_file, const int threads) {
    BOOST_LOG_TRIVIAL(info) << "Loading queries.";

    const auto features = wasabi::mass::load_featurexml(feature_file);

    const auto points = wasabi::mass::load_mzml(mzml_file);

    auto feature_specs = wasabi::mass::generate_compact_feature_entries(
        features, points, bin_mz_k, bin_width_k, threads
    );

    BOOST_LOG_TRIVIAL(info) << "Loaded " + std::to_string(feature_specs.size()) + " feature entries from "
        << features.size() << " features and " << points.size() << " points.";

    return feature_specs;
}

wasabi::db::CompactDistribution<wasabi::mz_t, wasabi::intensity_t> to_distribution2(const wasabi::mass::CompactFeatureSpec& f)
{
    return { f.size, f.values.get(), f.probabilities.get() };
}

template<typename T>
std::vector<size_t> explode(const std::vector<size_t>& indexes, const PeptideEntries& entries,
    const std::vector<std::tuple<T, size_t, size_t>>& mapping, const double rt, const double rt_threshold) {
    std::vector<size_t> ret;
    for (auto idx : indexes) {
        const size_t start = std::get<1>(mapping[idx]);
        const size_t n = std::get<2>(mapping[idx]);
        for (size_t i = start; i < start + n; ++i) {
            if (std::abs(entries[i].rt - rt) < rt_threshold) {
                ret.push_back(i);
            }
        }
    }
    return ret;
}

std::pair<std::vector<double>, std::vector<size_t>> split_items(const std::vector<std::pair<double, size_t>> &items) {
    std::vector<double> values(items.size());
    std::vector<size_t> indices(items.size());

    for (size_t i = 0; i < values.size(); ++i) {
        values[i] = items[i].first;
        indices[i] = items[i].second;
    }

    return {std::move(values), std::move(indices)};
}

search_result_t answer_queries(const database_t& db,
                               const PeptideEntries& entries,
                               const queries_t& queries,
                               const compressed_t& compressed_entries,
                               const double rt_threshold,
                               const int threads)
{
    BOOST_LOG_TRIVIAL(info) << "Searching the database for answers.";

    std::vector<std::vector<size_t>> best_indices(queries.size());
    std::vector<std::vector<double>> best_dists(queries.size());

#pragma omp parallel num_threads(threads) default(none) shared(queries, db, entries, compressed_entries, rt_threshold, best_indices, best_dists)
    {
#pragma omp for schedule(static)
        for (size_t i = 0; i < queries.size(); ++i) {
            const auto& f = queries[i];
            auto query_dist = to_distribution2(f);
            auto search_result = db.search(query_dist);
            const auto& [distvals, idxs] = split_items(search_result);

            // Explode indices AND align distances
            auto& exploded_indices =best_indices[i];
            auto& exploded_dists = best_dists[i];
            for (size_t k = 0; k < idxs.size(); ++k) {
                const size_t start = std::get<1>(compressed_entries[idxs[k]]);
                const size_t n     = std::get<2>(compressed_entries[idxs[k]]);
                for (size_t j = start; j < start + n; ++j) {
                    if (std::abs(entries[j].rt - f.rt) < rt_threshold) {
                        exploded_indices.push_back(j);
                        exploded_dists.push_back(distvals[k]);  // align distance
                    }
                }
            }
        }
    }

    return {std::move(best_indices), std::move(best_dists)};
}

void write_better_feature_annotations_to_csv(
    const std::string& filename,
    const std::vector<wasabi::mass::CompactFeatureSpec>& feature_specs,
    const std::vector<std::vector<size_t>>& best_indices,
    const std::vector<size_t>& annotations,
    const std::vector<wasabi::FastaRecord>& fasta_records,
    const std::vector<wasabi::iso::PeptideEntry>& entries,
    const std::vector<std::vector<double>>& dists)
{
    std::ofstream csv(filename);
    if (!csv)
        throw std::runtime_error("Cannot open CSV for writing: " + filename);

    csv << std::fixed
        << std::setprecision(std::numeric_limits<double>::max_digits10);

    // Header
    csv << "FEATURE_ID,SEQUENCE,ADDUCT,CHARGE,MODIFICATIONS,"
           "TOTAL_COMPOSITION,IS_DECOY,DIST,RT_pred\n";

    if (feature_specs.size() != best_indices.size())
        throw std::runtime_error("feature_specs and best_indices size mismatch");

    if (feature_specs.size() != dists.size())
        throw std::runtime_error("feature_specs and dists size mismatch");

    for (size_t i = 0; i < feature_specs.size(); ++i) {
        const auto& feature = feature_specs[i];

        for (size_t k = 0; k < best_indices[i].size(); ++k) {
            size_t annot_idx = best_indices[i][k];
            double dist = dists[i][k];

            if (annot_idx >= annotations.size()) {
                throw std::runtime_error(
                    "Annotation index out of bounds at feature_id " +
                    std::to_string(i));
            }

            const auto& record = fasta_records[annotations[annot_idx]];
            const auto& entry  = entries[annot_idx];

            // FEATURE_ID, SEQUENCE, ADDUCT, CHARGE
            csv << feature.id << ",";
            csv << record.sequence << ",";
            csv << *entry.adduct_name << ",";
            csv << std::to_string(entry.charge) << ",\"";

            bool first = true;

            if (entry.peptide) {
                // Fixed PTMs
                for (uint8_t j = 0; j < entry.peptide->residue_mods_size; j++) {
                    const auto& [count, ptmIdx] = entry.peptide->residue_mods_ptr[j];
                    const auto& ptm = wasabi::ptm::common_ptms_k[ptmIdx];
                    if (ptm.type != wasabi::ptm::PtmType::Fixed) continue;

                    if (!first) csv << ";";
                    csv << ptm.name;

                    size_t n = std::count(
                        record.sequence.begin(),
                        record.sequence.end(),
                        ptm.targetResidue);

                    if (n > 1) csv << "x" << n;
                    first = false;
                }

                // Variable PTMs
                for (uint8_t j = 0; j < entry.peptide->residue_mods_size; j++) {
                    const auto& [count, ptmIdx] = entry.peptide->residue_mods_ptr[j];
                    const auto& ptm = wasabi::ptm::common_ptms_k[ptmIdx];
                    if (ptm.type != wasabi::ptm::PtmType::Variable) continue;

                    if (!first) csv << ";";
                    csv << ptm.name;
                    if (count > 1) csv << "x" << count;
                    first = false;
                }

                if (entry.peptide->has_c_term()) {
                    if (!first) csv << ";";
                    csv << entry.peptide->get_c_term().name << "@C-term";
                    first = false;
                }

                if (entry.peptide->has_n_term()) {
                    if (!first) csv << ";";
                    csv << entry.peptide->get_n_term().name << "@N-term";
                }
            }

            csv << "\",";
            csv << entry.total_composition.to_string() << ",";
            csv << (entry.is_decoy ? "true" : "false") << ",";
            csv << dist << ",";
            csv << entry.rt;
            csv << "\n";
        }
    }
}

void save_annotations(const wasabi::FastaRecords& records,
                      const PeptideEntries& entries,
                      const std::string& output_file,
                      const queries_t& queries,
                      const search_result_t& search_result)
{
    BOOST_LOG_TRIVIAL(info) << "Saving annotations to \"" << output_file << "\"";

    const auto [indices, dists] = search_result;

    const auto better_annots = generate_better_annotations(records, entries);

    write_better_feature_annotations_to_csv(output_file,
                                            queries,
                                            indices,
                                            better_annots,
                                            records,
                                            entries,
                                            dists);
}
