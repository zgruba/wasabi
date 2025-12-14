#pragma once

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdexcept>
#include <cerrno>
#include <vector>
#include <fcntl.h>
#include <cstring>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <iomanip>
#include <queue>
#include <iostream>

namespace wasabi::db {

// ===============================================================
//                     TYPE DEFINITIONS
// ===============================================================

/** @brief Identifier for database entries. */
using entry_id_t = uint64_t;

/** @brief Number of elements in compact distribution. */
using distribution_size_t = uint16_t;

/** @brief Byte offset inside mmap'ed file. */
using offset_t = uint64_t;

template<typename K, typename V, typename = std::enable_if_t<std::is_arithmetic_v<K>>>
class MinThresholdKVPQ {
public:
    explicit MinThresholdKVPQ(K threshold) : threshold_(std::abs(threshold)) {}

    // push() assumes push_first() was called first.
    // Otherwise, min_key_ == max(), and this would reject initial items.
    void push(const K key, const V &value) {
        if (key <= min_key_ + threshold_) {
            pq_.push({key, value});
            if (key < min_key_) {
                min_key_ = key;
                remove_outliers();
            }
        }
    }

    // push_first() intentionally does NOT remove outliers.
    // Used for initial seeding or bulk ingest where filtering may be deferred.
    void push_first(const K key, const V &value) {
        min_key_ = key;
        pq_.push({key, value});
    }

    std::vector<V> extract_values() {
        std::vector<V> res;
        while (!pq_.empty()) {
            res.push_back(pq_.top().second);
            pq_.pop();
        }
        min_key_ = std::numeric_limits<K>::max();
        return res;
    }

    // New: extract only keys in priority order
    std::vector<K> extract_keys() {
        std::vector<K> res;
        while (!pq_.empty()) {
            res.push_back(pq_.top().first);
            pq_.pop();
        }
        min_key_ = std::numeric_limits<K>::max();
        return res;
    }

    std::vector<std::pair<K,V>> extract_items() {
        std::vector<std::pair<K,V>> res;
        while (!pq_.empty()) {
            res.push_back(pq_.top());
            pq_.pop();
        }
        min_key_ = std::numeric_limits<K>::max();
        return res;
    }

    K min_key() const { return min_key_; }

private:
    const K threshold_;
    K min_key_ = std::numeric_limits<K>::max();
    std::priority_queue<std::pair<K, V>> pq_{};

    void remove_outliers() {
        // assume !pq.empty() and threshold_ > 0
        while (std::abs(pq_.top().first - min_key_) > threshold_) {
            pq_.pop();
        }
    }
};

/**
 * @brief One element of a probability distribution.
 *
 * @tparam ValueT  Value type
 * @tparam ProbT   Probability type
 */
template <typename ValueT, typename ProbT>
struct DistributionElement {
    ValueT value; ///< Data point location.
    ProbT probability; ///< Probability mass at this location.
};

/**
 * @brief A dynamic in-RAM distribution stored in `std::vector`.
 */
template <typename ValueT, typename ProbT>
using Distribution = std::vector<DistributionElement<ValueT, ProbT>>;

/**
 * @brief Identifier type for an element inside a dynamic distribution.
 */
template <typename ValueT, typename ProbT>
using DistributionId = Distribution<ValueT, ProbT>::size_type;


/**
 * @brief A compact, pointer-based distribution used for mmap'ed file access.
 *
 * This avoids allocations and copying by directly pointing into the mmap'ed region.
 */
template <typename ValueT, typename ProbT>
struct CompactDistribution {
    using size_type = distribution_size_t;
    distribution_size_t size; ///< Number of points.
    ValueT* values;           ///< Pointer to sorted values.
    ProbT* probabilities;     ///< Pointer to corresponding probabilities.
};

/** @brief Identifier for CompactDistribution elements. */
template <typename ValueT, typename ProbT>
using CompactDistributionId = CompactDistribution<ValueT, ProbT>::size_type;

// ===============================================================
//                   WASSERSTEIN DISTANCE
// ===============================================================
/**
 * @brief Computes q-Wasserstein distance between two compact distributions.
 *
 * Both distributions must be sorted by value in ascending order.
 *
 * @tparam ValueT Numeric type of values.
 * @tparam ProbT  Numeric type of probabilities.
 * @param a First distribution.
 * @param b Second distribution.
 * @param q Exponent parameter (typically 1).
 * @return The q-Wasserstein distance.
 */
template <typename ValueT, typename ProbT>
__always_inline ValueT wasserstein_distance(const CompactDistribution<ValueT, ProbT>& a,
                            const CompactDistribution<ValueT, ProbT>& b,
                            const int q = 1) {
    distribution_size_t i = 0, j = 0;
    ProbT cdf_a = 0.0, cdf_b = 0.0;
    ValueT prev_x = 0.0, curr_x = 0.0, accumulator = 0.0;

    // Assume sizes are non 0;
    if (a.values[0] <= b.values[0]) {
        prev_x = a.values[0];
    } else {
        prev_x = b.values[0];
    }

    while (i < a.size && a.values[i] == prev_x) {
        cdf_a += a.probabilities[i];
        ++i;
    }

    while (j < b.size && b.values[j] == prev_x) {
        cdf_b += b.probabilities[j];
        ++j;
    }

    while (i < a.size || j < b.size) {
        if (i == a.size) {
            curr_x = b.values[j];
        } else if (j == b.size) {
            curr_x = a.values[i];
        } else {
            curr_x = std::min(a.values[i], b.values[j]);
        }

        auto dx = curr_x - prev_x;
        auto diff = std::abs(cdf_a - cdf_b);
        accumulator += std::pow(diff, q) * dx;

        while (i < a.size && a.values[i] == curr_x) {
            cdf_a += a.probabilities[i];
            ++i;
        }

        while (j < b.size && b.values[j] == curr_x) {
            cdf_b += b.probabilities[j];
            ++j;
        }

        prev_x = curr_x;
    }

    return std::pow(accumulator, 1.0 / q);
}

/**
 * @brief Computes q-Wasserstein distance between two vector-based distributions.
 *
 * The input vectors must be sorted by value.
 *
 * @tparam ValueT Numeric type of values.
 * @tparam ProbT  Numeric type of probabilities.
 * @param a First distribution.
 * @param b Second distribution.
 * @param q Exponent parameter (typically 1).
 */
template <typename ValueT, typename ProbT>
__always_inline ValueT wasserstein_distance(const std::vector<DistributionElement<ValueT, ProbT>>& a,
                                            const std::vector<DistributionElement<ValueT, ProbT>>& b,
                                            const int q = 1) {
    // Make common axis for mz
    std::vector<double> mz_points;
    mz_points.reserve(a.size() + b.size());
    for (auto&[mz, probability] : a) {
        mz_points.push_back(mz);
    }
    for (auto&[mz, probability] : b) {
        mz_points.push_back(mz);
    }
    std::ranges::sort(mz_points);
    mz_points.erase(std::ranges::unique(mz_points).begin(), mz_points.end());

    // Calculate cdf
    std::vector<double> cdf_a(mz_points.size(), 0.0), cdf_b(mz_points.size(), 0.0);
    double running_a = 0.0, running_b = 0.0;
    size_t idx_a = 0, idx_b = 0;

    for (size_t i = 0; i < mz_points.size(); ++i) {
        const double mz = mz_points[i];
        while (idx_a < a.size() && a[idx_a].value <= mz)
            running_a += a[idx_a++].probability;
        while (idx_b < b.size() && b[idx_b].value <= mz)
            running_b += b[idx_b++].probability;
        cdf_a[i] = running_a;
        cdf_b[i] = running_b;
    }

    // Wasserstein explicit formula
    double accum = 0.0;
    for (size_t i = 1; i < mz_points.size(); ++i) {
        const double dx = mz_points[i] - mz_points[i - 1];
        const double diff = std::abs(cdf_a[i - 1] - cdf_b[i - 1]);
        accum += std::pow(diff, q) * dx;
    }

    return std::pow(accum, 1.0 / q);
}

// ===============================================================
//                  MEMORY-MAPPED FILE DATABASE
// ===============================================================

/**
 * @brief A database storing distributions inside a memory-mapped file.
 *
 * Sorted by expected value and searching by Wasserstein distance.
 */
template <typename ValueT, typename ProbT>
class FileQWassersteinDatabase {
public:
    /**
     * @brief Header metadata for each mmap'ed entry.
     */
    struct HeaderElement {
        entry_id_t id;               ///< Unique identifier.
        ValueT expected;             ///< Precomputed expected value.
        distribution_size_t size;    ///< Number of distribution points.
        offset_t offset;             ///< Offset to start of {values, probabilities}.
    };

    /**
     * @brief Open and memory-map database file.
     * @param filename Path to binary database file.
     * @throws std::runtime_error If file cannot be opened or mmap'ed.
     */
    explicit FileQWassersteinDatabase(const std::string& filename);
    /**
     * @brief Unmap file and release resources.
     */
    ~FileQWassersteinDatabase();
    /**
     * @brief Search using vector-based query.
     */
    [[nodiscard]] DistributionId<ValueT, ProbT> search(const Distribution<ValueT, ProbT>& query) const;
    /**
    * @brief Search using compact pointer-based query.
    *
    * Much faster because no allocations happen.
    */
    [[nodiscard]] CompactDistributionId<ValueT, ProbT> searchCompact(const CompactDistribution<ValueT, ProbT>& query) const;

private:
    struct Entry {
        entry_id_t id;
        ValueT expected;
        distribution_size_t size;
        ValueT* mzs;
        ProbT* probs;

        explicit Entry(const HeaderElement& e, void* start) : id(e.id), expected(e.expected), size(e.size) {
            mzs = reinterpret_cast<ValueT*>(reinterpret_cast<offset_t>(start) + e.offset);
            probs = static_cast<ProbT*>(mzs + e.size);
        }
    };
    uint64_t num_entries_;         
    std::vector<Entry> entries_;   

    const std::string& filename_;  
    int fd_ = -1;                  
    void* data_ = nullptr;         
    size_t size_ = 0;              
    
    static CompactDistribution<ValueT, ProbT> to_compact_distribution(const Entry& e) noexcept {
        return { e.size, e.mzs, e.probs };
    }
};

/**
 * @brief Open and memory-map a Wasserstein database file.
 *
 * This constructor opens the specified binary file, reads the header metadata,
 * and memory-maps the file into RAM. Each entry
 * in the file is wrapped as an internal Entry object.
 *
 * @param filename Path to the binary database file.
 * @throws std::runtime_error if the file cannot be opened, stat-ed, or mmap-ed.
 */
template <typename ValueT, typename ProbT>
FileQWassersteinDatabase<ValueT, ProbT>::FileQWassersteinDatabase(const std::string& filename) : filename_(filename) {
    fd_ = open(filename.c_str(), O_RDONLY);
    if (fd_ == -1) {
        throw std::runtime_error("open() failed for " + filename + ": " + std::strerror(errno));
    }

    struct stat st{};
    if (fstat(fd_, &st) == -1) {
        close(fd_);
        throw std::runtime_error("fstat() failed for " + filename + ": " + std::strerror(errno));
    }
    size_ = static_cast<size_t>(st.st_size);

    data_ = ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (data_ == MAP_FAILED) {
        close(fd_);
        throw std::runtime_error("mmap() failed for " + filename + ": " + std::strerror(errno));
    }

    close(fd_);
    fd_ = -1;

    const auto num_entries_ptr = static_cast<entry_id_t*>(data_);
    auto headers = static_cast<HeaderElement*>(static_cast<void*>(num_entries_ptr + 1));
    num_entries_ = *num_entries_ptr;
    entries_.reserve(num_entries_);
    for (int i = 0; i< num_entries_; ++i) {
        entries_.emplace_back(headers[i], data_);
    }
}

/**
 * @brief Unmap the database file from memory and release resources.
 *
 * This destructor calls `munmap` on the memory-mapped region to free resources.
 */
template <typename ValueT, typename ProbT>
FileQWassersteinDatabase<ValueT, ProbT>::~FileQWassersteinDatabase() {
    munmap(data_, size_);
}

/**
 * @brief Convert a raw pointer-based distribution into a vector-based one.
 *
 * This helper constructs a dynamic `Distribution` (AoS layout) by copying
 * `n` pairs of {value, probability} from raw pointer arrays.
 *
 * @warning This performs an allocation and `n` element copies. It is intended
 *          only for convenience and is significantly slower than using
 *          CompactDistribution directly.
 *
 * @param n Number of elements in the distribution.
 * @param v Pointer to contiguous values  (length ≥ n).
 * @param p Pointer to contiguous probabilities (length ≥ n).
 * @return A newly allocated `Distribution<ValueT, ProbT>` containing the data.
 */
template <typename ValueT, typename ProbT>
Distribution<ValueT, ProbT> p2d(const uint64_t n, ValueT* v, ProbT* p) {
    Distribution<ValueT, ProbT> distribution;
    for (uint64_t i = 0; i < n; ++i) {
        DistributionElement<ValueT, ProbT> x = {v[i], p[i]};
        distribution.push_back(x);
    }
    return distribution;
}

/**
 * @brief Search for the nearest distribution using a vector-based query.
 *
 * Computes the expected value of the query distribution, performs a binary search
 * to find the closest expected value in the database, and then refines the result
 * using the q-Wasserstein distance. Checks neighboring entries to ensure the
 * closest match is found.
 *
 * @param query The query distribution stored as a vector of DistributionElement.
 * @return The identifier of the nearest distribution in the memory-mapped database.
 */
template <typename ValueT, typename ProbT>
DistributionId<ValueT, ProbT>
FileQWassersteinDatabase<ValueT, ProbT>::search(const Distribution<ValueT, ProbT>& query) const {
    double min_diff = std::numeric_limits<double>::infinity();
    Distribution<ValueT, ProbT> best_dist;
    size_t best_index = 0;

    // Compute expected value
    ValueT expected_value{};
    ProbT total_prob{};

    for (const auto& elem : query) {
        expected_value += elem.value * elem.probability;
        total_prob += elem.probability;
    }

    if (total_prob != static_cast<ProbT>(0)) {
        expected_value /= total_prob;
    }

    // --- O(log n) binary search for closest expected value ---
    auto it = std::upper_bound(entries_.begin(), entries_.end(), expected_value,
        [](double value, const auto& e) {
            return value < e.expected;
        });

    if (it != entries_.begin()) {
        auto prev_it = std::prev(it);
        if (const double diff = std::abs(prev_it->expected - expected_value); diff < min_diff) {
            min_diff = diff;
            best_dist = p2d(prev_it->size, prev_it->mzs, prev_it->probs);
            best_index = std::distance(entries_.begin(), prev_it);
        }
    }
    if (it != entries_.end()) {
        double diff = std::abs(it->expected - expected_value);
        if (diff < min_diff) {
            min_diff = diff;
            best_dist = p2d(it->size, it->mzs, it->probs);
            best_index = std::distance(entries_.begin(), it);
        }
    }
    // ----------------------------------------------------------

    // Calculate Wasserstein distance
    double wmin = wasserstein_distance(query, best_dist);
    const size_t idx = best_index;

    // Check neighbours to the left
    for (int i = static_cast<int>(idx) - 1; i >= 0; --i) {
        if (const double Wlb = std::abs(expected_value - entries_[i].expected); Wlb > wmin) break;
        if (const double W = wasserstein_distance(query, p2d(entries_[i].size, entries_[i].mzs, entries_[i].probs)); W < wmin) {
            wmin = W;
            best_index = i;
        }
    }

    // Check neighbours to the right
    for (size_t i = idx + 1; i < entries_.size(); ++i) {
        if (const double Wlb = std::abs(expected_value - entries_[i].expected); Wlb > wmin) break;
        if (const double W = wasserstein_distance(query, p2d(entries_[i].size, entries_[i].mzs, entries_[i].probs)); W < wmin) {
            wmin = W;
            best_index = i;
        }
    }

    return entries_[best_index].id;
}

/**
 * @brief Search for the nearest distribution using a compact pointer-based query.
 *
 * Similar to `search`, but operates on a CompactDistribution to avoid vector
 * allocations and copies. Computes the expected value,
 * performs binary search for closest expected value, and refines using
 * the q-Wasserstein distance, checking neighbours for better matches.
 *
 * @param query The query distribution stored as a CompactDistribution.
 * @return The identifier of the nearest distribution in the memory-mapped database.
 */
template <typename ValueT, typename ProbT>
CompactDistributionId<ValueT, ProbT>
FileQWassersteinDatabase<ValueT, ProbT>::searchCompact(const CompactDistribution<ValueT, ProbT>& query) const
{
    double min_diff = std::numeric_limits<double>::infinity();
    size_t best_index = 0;

    // --- Compute expected value of the query distribution ---
    ValueT expected_value = 0;
    ProbT total_prob = 0;

    for (distribution_size_t i = 0; i < query.size; ++i) {
        expected_value += query.values[i] * query.probabilities[i];
        total_prob += query.probabilities[i];
    }

    if (total_prob != static_cast<ProbT>(0)) {
        expected_value /= total_prob;
    }

    // --- Binary search for closest expected value ---
    auto it = std::upper_bound(entries_.begin(), entries_.end(), expected_value,
        [](double value, const auto& e) {
            return value < e.expected;
        });

    if (it != entries_.begin()) {
        auto prev_it = std::prev(it);
        const double diff = std::abs(prev_it->expected - expected_value);
        if (diff < min_diff) {
            min_diff = diff;
            best_index = std::distance(entries_.begin(), prev_it);
        }
    }
    if (it != entries_.end()) {
        const double diff = std::abs(it->expected - expected_value);
        if (diff < min_diff) {
            min_diff = diff;
            best_index = std::distance(entries_.begin(), it);
        }
    }

    // --- Refine using Wasserstein distance around the best match ---
    double wmin = wasserstein_distance(query, to_compact_distribution(entries_[best_index]), 1);
    const size_t idx = best_index;

    // Search left side
    for (int i = static_cast<int>(idx) - 1; i >= 0; --i) {
        if (const double Wlb = std::abs(expected_value - entries_[i].expected); Wlb > wmin) break;

        if (const double W = wasserstein_distance(query, to_compact_distribution(entries_[i]), 1); W < wmin) {
            wmin = W;
            best_index = i;
        }
    }

    // Search right side
    for (size_t i = idx + 1; i < entries_.size(); ++i) {
        if (const double Wlb = std::abs(expected_value - entries_[i].expected); Wlb > wmin) break;

        if (const double W = wasserstein_distance(query, to_compact_distribution(entries_[i]), 1); W < wmin) {
            wmin = W;
            best_index = i;
        }
    }

    return entries_[best_index].id;
}


// ===============================================================
//                IN-RAM SORTED WASSERSTEIN DATABASE
// ===============================================================

/**
 * @brief In-RAM database with entries sorted by expected value.
 */
template <typename ValueT, typename ProbT>
class QWassersteinDatabase {
    struct Entry {
        ValueT expected_value;                         ///< Expected value of distribution.
        Distribution<ValueT, ProbT> distribution;      ///< Full distribution.
        DistributionId<ValueT, ProbT> distribution_id; ///< Index.

        long long getSize() {
            auto s1 = distribution.capacity() * (sizeof(decltype(distribution.at(0))));
            return sizeof(expected_value) + s1 + sizeof(distribution_id);
        }
    };
    std::vector<Entry> entries; ///< Sorted by expected_value.
public:
    /**
     * @brief Build database from list of distributions.
     */
    explicit QWassersteinDatabase(const std::vector<Distribution<ValueT, ProbT>> &distributions);

    /**
     * @brief Nearest-neighbor search using Wasserstein distance.
     */
    [[nodiscard]] DistributionId<ValueT, ProbT> search(const Distribution<ValueT, ProbT>& query) const;
};

/**
 * @brief Construct an in-RAM Wasserstein database from a list of distributions.
 *
 * Computes the expected value for each distribution, stores it along with the
 * full distribution and its index, and sorts all entries by expected value.
 *
 * @param distributions Vector of distributions to store in the database.
 */
template <typename ValueT, typename ProbT>
QWassersteinDatabase<ValueT, ProbT>::QWassersteinDatabase(const std::vector<Distribution<ValueT, ProbT>> &distributions) {
    this->entries.clear();
    this->entries.reserve(distributions.size());

    // Create entries from input distributions
    for (size_t i = 0; i < distributions.size(); ++i) {
        const auto& dist = distributions[i];

        // Compute expected value
        ValueT expected_value{};
        ProbT total_prob{};

        for (const auto& elem : dist) {
            expected_value += elem.value * elem.probability;
            total_prob += elem.probability;
        }

        if (total_prob != static_cast<ProbT>(0)) {
            expected_value /= total_prob;
        }

        // Construct entry
        Entry entry;
        entry.expected_value = expected_value;
        entry.distribution = dist;
        entry.distribution_id = static_cast<DistributionId<ValueT, ProbT>>(i);

        this->entries.push_back(std::move(entry));
    }
    // Ensure entries are sorted by expected_mass
    std::sort(this->entries.begin(), this->entries.end(),
              [](const Entry& a, const Entry& b) {
                  return a.expected_value < b.expected_value;
              });

    long long v = 0;
    for (size_t i = 1; i < this->entries.size(); ++i) {
        v += entries[i].getSize();
    }
    std::cout << "sizeof(Entry) = " << v << std::endl;
}

/**
 * @brief Search for the nearest distribution using the q-Wasserstein distance.
 *
 * Computes the expected value of the query distribution, performs a binary search
 * to find the closest expected value in the database, and then refines the search
 * by comparing Wasserstein distances with neighboring entries.
 *
 * @param query The query distribution stored as a vector of DistributionElement.
 * @return The identifier of the nearest distribution in the database.
 */
template <typename ValueT, typename ProbT>
DistributionId<ValueT, ProbT>QWassersteinDatabase<ValueT, ProbT>::search(const Distribution<ValueT, ProbT>& query) const {
    double min_diff = std::numeric_limits<double>::infinity();
    Distribution<ValueT, ProbT> best_dist;
    size_t best_index = 0;

    // --- Compute expected value of the query distribution ---
    ValueT expected_value{};
    ProbT total_prob{};

    for (const auto& elem : query) {
        expected_value += elem.value * elem.probability;
        total_prob += elem.probability;
    }

    if (total_prob != static_cast<ProbT>(0)) {
        expected_value /= total_prob;
    }

    // --- O(log n) binary search for closest expected_value ---
    auto it = std::upper_bound(entries.begin(), entries.end(), expected_value,
        [](double value, const auto& e) {
            return value < e.expected_value;
        });

    if (it != entries.begin()) {
        auto prev_it = std::prev(it);
        if (const double diff = std::abs(prev_it->expected_value - expected_value); diff < min_diff) {
            min_diff = diff;
            best_dist = prev_it->distribution;
            best_index = std::distance(entries.begin(), prev_it);
        }
    }

    if (it != entries.end()) {
        if (const double diff = std::abs(it->expected_value - expected_value); diff < min_diff) {
            min_diff = diff;
            best_dist = it->distribution;
            best_index = std::distance(entries.begin(), it);
        }
    }
    // ---------------------------------------------------------

    // --- Refine using Wasserstein distance around the best match ---
    double wmin = wasserstein_distance(query, best_dist);
    const size_t idx = best_index;

    // Search left side (smaller expected values)
    for (int i = static_cast<int>(idx) - 1; i >= 0; --i) {
        if (const double Wlb = std::abs(expected_value - entries[i].expected_value); Wlb > wmin) break;
        if (const double W = wasserstein_distance(query, entries[i].distribution); W <= wmin) {
            wmin = W;
            best_index = i;
        }
    }

    // Search right side (larger expected values)
    for (size_t i = idx + 1; i < entries.size(); ++i) {
        if (const double Wlb = std::abs(expected_value - entries[i].expected_value); Wlb > wmin) break;
        if (const double W = wasserstein_distance(query, entries[i].distribution); W <= wmin) {
            wmin = W;
            best_index = i;
        }
    }

    return entries[best_index].distribution_id;
}


// ===============================================================
//                 NAIVE O(N) WASSERSTEIN DATABASE
// ===============================================================

/**
 * @brief A slow but simple baseline database (linear scan).
 */
template <typename ValueT, typename ProbT>
class NaiveWassersteinDatabase {
    struct Entry {
        Distribution<ValueT, ProbT> distribution;
        DistributionId<ValueT, ProbT> distribution_id;
    };
    std::vector<Entry> entries; // kept sorted by expected value
public:
    /**
     * @brief Construct naive database.
     */
    explicit NaiveWassersteinDatabase(const std::vector<Distribution<ValueT, ProbT>> &distributions);
    /**
     * @brief Linear scan for nearest distribution.
     */
    [[nodiscard]] DistributionId<ValueT, ProbT> search(const Distribution<ValueT, ProbT>& query) const;
};

/**
 * @brief Construct a naive Wasserstein database from a list of distributions.
 *
 * Computes the expected value for each distribution, stores it along with the
 * full distribution and its index. This database uses a linear scan for searching.
 *
 * @param distributions Vector of distributions to store in the database.
 */
template <typename ValueT, typename ProbT>
NaiveWassersteinDatabase<ValueT, ProbT>::NaiveWassersteinDatabase(const std::vector<Distribution<ValueT, ProbT>> &distributions) {
    this->entries.clear();
    this->entries.reserve(distributions.size());

    // Create entries from input distributions
    for (size_t i = 0; i < distributions.size(); ++i) {
        const auto& dist = distributions[i];

        // Compute expected value
        ValueT expected_value{};
        ProbT total_prob{};

        for (const auto& elem : dist) {
            expected_value += elem.value * elem.probability;
            total_prob += elem.probability;
        }

        if (total_prob != static_cast<ProbT>(0)) {
            expected_value /= total_prob;
        }

        // Construct entry
        Entry entry;
        entry.expected_value = expected_value;
        entry.distribution = dist;
        entry.distribution_id = static_cast<DistributionId<ValueT, ProbT>>(i);

        this->entries.push_back(std::move(entry));
    }
}

/**
 * @brief Linear search for the nearest distribution using Wasserstein distance.
 *
 * Iterates over all entries in the database and computes the Wasserstein distance
 * to the query distribution, returning the identifier of the closest match.
 *
 * @param query The query distribution stored as a vector of DistributionElement.
 * @return The identifier of the nearest distribution in the database.
 * @throws std::runtime_error if the database is empty.
 */
template <typename ValueT, typename ProbT>
DistributionId<ValueT, ProbT> NaiveWassersteinDatabase<ValueT, ProbT>::search(const Distribution<ValueT, ProbT>& query) const {
    if (entries.empty()) {
        throw std::runtime_error("Database is empty");
    }

    double wmin = std::numeric_limits<double>::max();
    Entry best_entry = entries[0];

    for (const auto& e : entries) {
        if (const double W = wasserstein_distance(query, e.distribution); W < wmin) {
            wmin = W;
            best_entry = e;
        }
    }
    return best_entry.distribution_id;
}

// ===============================================================
//               STREAMING / ON-DEMAND DATABASE
// ===============================================================

/**
 * @brief Streaming database: distributions are generated on demand.
 *
 * This avoids storing all distributions in memory and uses only
 * precomputed expected values + a distribution generator.
 */
template <typename ValueT, typename ProbT>
class StreamQWassersteinDatabase {
public:
    /**
     * @brief Function generating distribution #idx.
     */
    using DistGenerator = std::function<CompactDistribution<ValueT, ProbT>(size_t idx)>;

    /**
     * @brief Construct streaming database.
     *
     * @param evs Vector of (expected value, index)
     * @param generator Callback generating distribution data
     * @param queue_threshold Threshold for storing elements in the queue
     * @param search_threshold Threshold for the search stop
     */
    StreamQWassersteinDatabase(
        const std::vector<std::pair<ValueT, size_t>>& evs,
        const DistGenerator& generator,
        const double queue_threshold,
        const double search_threshold
   ) : evs_(evs), dist_generator_(generator), queue_threshold_(queue_threshold), search_threshold_(search_threshold) {}

    /**
     * @brief Search for nearest neighbours.
     *
     * Returns:
     *   - a vector of pairs (Wdist, index)
     *   for nearest neighbours in selected thresholds.
     */
    std::vector<std::pair<ValueT,size_t>> search(const CompactDistribution<ValueT, ProbT>& query) const {
        // --- Compute EV of query ---
        ValueT query_ev = 0;
        ProbT total_prob = 0;
        for (size_t k = 0; k < query.size; ++k) {
            query_ev += query.values[k] * query.probabilities[k];
            total_prob += query.probabilities[k];
        }
        if (total_prob != static_cast<ProbT>(0)) query_ev /= total_prob;
        // --- Binary search on EV ---
        auto bound_it = std::upper_bound(
           evs_.begin(), evs_.end(), query_ev,
                                  [](ValueT val, const std::pair<ValueT, size_t>& e){ return val < e.first; }
        );
        size_t best_idx = 0;
        ValueT min_diff = std::numeric_limits<ValueT>::infinity();
        decltype(evs_.begin()) best_evs_it = bound_it;

        if (bound_it != evs_.begin()) {
            auto prev = std::prev(bound_it);
            ValueT diff = std::abs(prev->first - query_ev);
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = prev->second;
                best_evs_it = prev;
            }
        }
        if (bound_it != evs_.end()) {
            ValueT diff = std::abs(bound_it->first - query_ev);
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = bound_it->second;
                best_evs_it = bound_it;
            }
        }

        // --- Initial Wasserstein ---
        ValueT wmin = wasserstein_distance(query, dist_generator_(best_idx));
        std::vector<size_t> best_indexes;
        best_indexes.push_back(best_idx);

        auto best_indexes_pq = MinThresholdKVPQ<ValueT, size_t>(queue_threshold_);
        best_indexes_pq.push_first(wmin, best_idx);

        // --- Left refinement ---
        for (auto it = std::make_reverse_iterator(best_evs_it); it != evs_.rend(); ++it) {
            if (std::abs(query_ev - it->first) > wmin + search_threshold_) {
                break;
            }

            ValueT W = wasserstein_distance(query, dist_generator_(it->second));
            best_indexes_pq.push(W, it->second);
            wmin = best_indexes_pq.min_key();
        }

       // --- Right refinement ---
       for (auto it = best_evs_it + 1; it != evs_.end(); ++it) {
           if (std::abs(query_ev - it->first) > wmin + search_threshold_) {
               break;
           }

           ValueT W = wasserstein_distance(query, dist_generator_(it->second));
           best_indexes_pq.push(W, it->second);
           wmin = best_indexes_pq.min_key();
       }

        return best_indexes_pq.extract_items();

   }
    /**
     * @brief Dump distance profile to CSV.
     */
    void printDistribution(const std::string& filename, const CompactDistribution<ValueT, ProbT>& query) const {
        std::ofstream csv(filename);
        if (!csv)
            throw std::runtime_error("Cannot open CSV for writing: " + filename);

        // Header
        csv << "EV,Wdist\n" << std::setprecision(15);
        for (size_t k = 0; k < evs_.size(); ++k) {
            auto ev = evs_[k].first;
            auto idx = evs_[k].second;
            ValueT W = wasserstein_distance(query, dist_generator_(idx));
            csv << ev << "," << W << "\n";
        }
   }

private:
    const std::vector<std::pair<ValueT, size_t>>& evs_;
    const DistGenerator dist_generator_;
    const ValueT queue_threshold_;
    const ValueT search_threshold_;
};

} // namespace wasabi::db
