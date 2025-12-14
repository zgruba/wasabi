#include <string>
#include <vector>
#include <iostream>
#include <boost/program_options.hpp>
#include "iso.hpp"
#include "parsing.hpp"

namespace po = boost::program_options;

constexpr int missed_cleavages_k = 2;
constexpr int sequence_min_len_k = 7;
constexpr int sequence_max_len_k = 51;
constexpr bool allow_semitryptic_k = false;
constexpr int max_variable_modifications_k = 3;

struct Config {
    std::string path_to_fasta;
    std::string output_filename;
    int threads = 1;
    bool with_decoys = false;
};

void print_usage_and_exit(const po::options_description& all, int exit_code);

Config parse_args(int argc, char* argv[]);

void write_modifications(std::ofstream& out, const wasabi::iso::ModifiedPeptide& peptide, const std::string& sequence);

void extract_seq_and_mods_to_csv(const std::vector<wasabi::iso::PeptideEntry>& entries, const std::string& output_filename);

int main(const int argc, char* argv[]) {
    const auto [path_to_fasta, output_filename, threads, with_decoys] = parse_args(argc, argv);

    const auto records = wasabi::fasta::read_fasta_with_trypsin(
        path_to_fasta, missed_cleavages_k,sequence_min_len_k,sequence_max_len_k,allow_semitryptic_k, with_decoys);

    const auto entries = wasabi::iso::generate_entries_from_fasta_records(
        records, {1}, max_variable_modifications_k, threads);

    extract_seq_and_mods_to_csv(entries, output_filename);
}


void print_usage_and_exit(const po::options_description& all, int exit_code) {
    std::cout << "Usage:\n"
              << "  prepare_csv_for_deeplc <fasta> <output.csv> [options]\n\n"
              << all << "\n";
    std::exit(exit_code);
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;

    // Named Options:
    po::options_description named_opts("Options");
    named_opts.add_options()
        ("help,h", "Show help message")
        ("threads,t", po::value<int>()->default_value(1),
            "Number of threads (default: 1)")
        ("with-decoys,d", po::bool_switch(),
            "Include decoys (toggle)");

    // Positional placeholders (only for help):
    po::options_description positional_opts("Positional arguments");
    positional_opts.add_options()
        ("path_to_fasta", po::value<std::string>(), "Path to FASTA file (required)")
        ("output_filename", po::value<std::string>(), "Output CSV filename (required)");

    po::options_description all;
    all.add(named_opts).add(positional_opts);

    po::positional_options_description pos;
    pos.add("path_to_fasta", 1);
    pos.add("output_filename", 1);

    po::variables_map vm;

    // Parse arguments:
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all)
                      .positional(pos)
                      .run(),
                  vm);
    } catch (...) {
        // Parsing failed (bad option, missing value etc.) — show help and exit(1)
        print_usage_and_exit(all, 1);
    }

    if (vm.contains("help")) {
        print_usage_and_exit(all, 0);
    }

    // Validate presence of required positional args:
    if (!vm.contains("path_to_fasta") || !vm.contains("output_filename")) {
        std::cerr << "Error: missing required positional arguments.\n\n";
        print_usage_and_exit(all, 1);
    }

    // Assign values:
    cfg.path_to_fasta = vm["path_to_fasta"].as<std::string>();
    cfg.output_filename = vm["output_filename"].as<std::string>();
    cfg.threads = vm["threads"].as<int>();
    cfg.with_decoys = vm["with-decoys"].as<bool>();

    // Basic validation:
    if (cfg.threads <= 0) {
        std::cerr << "Error: threads must be >= 1.\n\n";
        print_usage_and_exit(all, 1);
    }

    return cfg;
}

void write_modifications(std::ofstream& out,
                         const wasabi::iso::ModifiedPeptide& peptide,
                         const std::string& sequence) {
    out << "\"";
    bool first = true;

    // Fixed PTMs:
    for (uint8_t j = 0; j < peptide.residue_mods_size; j++) {
        const auto& [count, ptmIdx] = peptide.residue_mods_ptr[j];
        if (const auto& ptm = wasabi::ptm::common_ptms_k[ptmIdx]; ptm.type == wasabi::ptm::PtmType::Fixed) {
            if (first) {
                first = false;
            } else {
                out << ";";
            }
            out << ptm.name;
            if (const size_t n = std::ranges::count(sequence, ptm.targetResidue); n > 1) {
                out << "x" << std::to_string(n);
            }
        }
    }

    // Variable PTMs:
    for (uint8_t j = 0; j < peptide.residue_mods_size; j++) {
        const auto& [count, ptmIdx] = peptide.residue_mods_ptr[j];
        if (const auto& ptm = wasabi::ptm::common_ptms_k[ptmIdx]; ptm.type == wasabi::ptm::PtmType::Variable) {
            if (first) {
                first = false;
            } else {
                out << ";";
            }
            out << ptm.name;
            if (count > 1) {
                out << "x" << std::to_string(count);
            }
        }
    }

    // Terminal mods:
    if (peptide.has_c_term()) {
        if (first) {
            first = false;
        } else {
            out << ";";
        }
        out << peptide.get_c_term().name << "@C-term";
    }
    if (peptide.has_n_term()) {
        if (!first) {
            out << ";";
        }
        out << peptide.get_n_term().name << "@N-term";
    }
    out << "\"";
}

void extract_seq_and_mods_to_csv(const std::vector<wasabi::iso::PeptideEntry>& entries,
                                 const std::string& output_filename) {
    constexpr auto HEADER_LINE = "SEQUENCE,MODIFICATIONS\n";
    std::ofstream csv(output_filename);
    if (!csv) {
        throw std::runtime_error("Cannot open file for writing: " + output_filename);
    }

    csv << HEADER_LINE;
    for (const auto& entry : entries) {
        csv << *entry.peptide->sequence << ",";
        write_modifications(csv, *entry.peptide, *entry.peptide->sequence);
        csv << "\n";
    }
}
