#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <mutex>
#include <atomic>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <deque>
#include <queue>
#include <zlib.h>

#include "BooPHF.h"

using kmer_t = uint64_t;

class TaxonomyTree {
public:
    TaxonomyTree() = default;
    
    void load_from_nodes_dmp(const std::string& nodes_path) {        
        std::ifstream file(nodes_path);
        if (!file) {
            throw std::runtime_error("Cannot open nodes.dmp: " + nodes_path);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            // nodes.dmp format: taxid | parent_taxid | rank | ...
            std::istringstream iss(line);
            uint32_t taxid, parent_taxid;
            char sep;
            
            if (iss >> taxid >> sep >> parent_taxid) {
                parent_map_[taxid] = parent_taxid;
                children_map_[parent_taxid].push_back(taxid);
            }
        }
        
    }
    
    std::unordered_set<uint32_t> get_all_descendants(uint32_t taxid) const {
        std::unordered_set<uint32_t> descendants;
        std::queue<uint32_t> to_visit;
        to_visit.push(taxid);
        
        while (!to_visit.empty()) {
            uint32_t current = to_visit.front();
            to_visit.pop();
            
            descendants.insert(current);
            
            auto it = children_map_.find(current);
            if (it != children_map_.end()) {
                for (uint32_t child : it->second) {
                    if (descendants.find(child) == descendants.end()) {
                        to_visit.push(child);
                    }
                }
            }
        }
        
        return descendants;
    }
    
    // Check if any of the descendants (excluding self) are in the target set
    bool has_target_descendants(uint32_t taxid, const std::unordered_set<uint32_t>& targets) const {
        std::queue<uint32_t> to_visit;
        std::unordered_set<uint32_t> visited;
        
        // Start with children, not self
        auto it = children_map_.find(taxid);
        if (it != children_map_.end()) {
            for (uint32_t child : it->second) {
                to_visit.push(child);
            }
        }
        
        while (!to_visit.empty()) {
            uint32_t current = to_visit.front();
            to_visit.pop();
            
            if (visited.find(current) != visited.end()) continue;
            visited.insert(current);
            
            // Check if this descendant is a target
            if (targets.find(current) != targets.end()) {
                return true;
            }
            
            // Add children to queue
            auto child_it = children_map_.find(current);
            if (child_it != children_map_.end()) {
                for (uint32_t child : child_it->second) {
                    to_visit.push(child);
                }
            }
        }
        
        return false;
    }
    
    uint32_t get_parent(uint32_t taxid) const {
        auto it = parent_map_.find(taxid);
        return (it != parent_map_.end()) ? it->second : 0;
    }
    
private:
    std::unordered_map<uint32_t, uint32_t> parent_map_;  // taxid -> parent
    std::unordered_map<uint32_t, std::vector<uint32_t>> children_map_;  // taxid -> children
};

class TargetTaxIDManager {
public:
    TargetTaxIDManager(const std::unordered_set<uint32_t>& targets, const TaxonomyTree& tree) {
        std::cout << "\nAnalyzing target taxid hierarchy..." << std::endl;
        
        // Categorize targets
        for (uint32_t target : targets) {
            if (tree.has_target_descendants(target, targets)) {
                // This target has other targets as descendants - exact match only
                targets_with_target_descendants_.insert(target);
                exact_match_taxids_.insert(target);
                std::cout << "  TaxID " << target << ": has target descendants (exact match only)" << std::endl;
            } else {
                // This target has NO other targets as descendants - rollup all descendants
                targets_without_target_descendants_.insert(target);
                
                // Get all descendants and map them to this target
                auto descendants = tree.get_all_descendants(target);
                for (uint32_t desc : descendants) {
                    descendant_to_target_[desc] = target;
                }
                
                std::cout << "  TaxID " << target << ": no target descendants (rolling up " 
                          << descendants.size() << " descendant taxids)" << std::endl;
            }
        }
        
        std::cout << "\nSummary:" << std::endl;
        std::cout << "  Targets with target descendants (exact match): " 
                  << targets_with_target_descendants_.size() << std::endl;
        std::cout << "  Targets without target descendants (rollup): " 
                  << targets_without_target_descendants_.size() << std::endl;
        std::cout << "  Total descendant mappings: " << descendant_to_target_.size() << std::endl;
    }

    std::pair<bool, uint32_t> get_target_for_taxid(uint32_t taxid) const {
        // Case 1: Exact match for targets with target descendants
        if (exact_match_taxids_.find(taxid) != exact_match_taxids_.end()) {
            return {true, taxid};
        }
        
        // Case 2: Check if it's a descendant of a rollup target
        auto it = descendant_to_target_.find(taxid);
        if (it != descendant_to_target_.end()) {
            return {true, it->second};
        }
        
        // Not a target and not a descendant of any rollup target
        return {false, 0};
    }
    
    // Get all taxids that should be considered (for the extraction phase)
    std::unordered_set<uint32_t> get_all_relevant_taxids() const {
        std::unordered_set<uint32_t> all;
        
        // Add exact match targets
        for (uint32_t t : exact_match_taxids_) {
            all.insert(t);
        }
        
        // Add all descendants of rollup targets
        for (const auto& [desc, target] : descendant_to_target_) {
            all.insert(desc);
        }
        
        return all;
    }
    
private:
    std::unordered_set<uint32_t> targets_with_target_descendants_;
    std::unordered_set<uint32_t> targets_without_target_descendants_;
    std::unordered_set<uint32_t> exact_match_taxids_;
    std::unordered_map<uint32_t, uint32_t> descendant_to_target_;  // any descendant -> target taxid
};

class MinimizerScanner {
public:
    MinimizerScanner(int k, int l, uint64_t spaced_seed_mask, uint64_t toggle_mask)
        : k_(k), l_(l), spaced_seed_mask_(spaced_seed_mask), 
          toggle_mask_(toggle_mask), str_pos_(0), loaded_ch_(0)
    {
        if (l > 31) {
            throw std::invalid_argument("l must be <= 31");
        }
        
        lmer_mask_ = (1ULL << (l * 2)) - 1;
        toggle_mask_ &= lmer_mask_;
        
        // Initialize lookup table
        for (int i = 0; i < 256; i++) {
            lookup_table_[i] = UINT8_MAX;
        }
        
        // Set DNA codes (case-insensitive)
        set_lookup('A', 0x00); set_lookup('a', 0x00);
        set_lookup('C', 0x01); set_lookup('c', 0x01);
        set_lookup('G', 0x02); set_lookup('g', 0x02);
        set_lookup('T', 0x03); set_lookup('t', 0x03);
    }
    
    void load_sequence(const std::string& seq) {
        str_ = &seq;
        str_pos_ = 0;
        str_len_ = seq.length();
        queue_.clear();
        queue_pos_ = 0;
        loaded_ch_ = 0;
        lmer_ = 0;
        last_minimizer_ = ~0ULL;
    }
    
    uint64_t* next_minimizer() {
        if (str_pos_ >= str_len_) return nullptr;
        
        bool changed_minimizer = false;
        
        while (!changed_minimizer) {
            // Incorporate next character
            if (loaded_ch_ == l_) {
                loaded_ch_--;
            }
            
            while (loaded_ch_ < l_ && str_pos_ < str_len_) {
                loaded_ch_++;
                lmer_ <<= 2;
                
                uint8_t lookup_code = lookup_table_[(unsigned char)(*str_)[str_pos_++]];
                
                if (lookup_code == UINT8_MAX) {
                    // Ambiguous character - reset
                    queue_.clear();
                    queue_pos_ = 0;
                    lmer_ = 0;
                    loaded_ch_ = 0;
                } else {
                    lmer_ |= lookup_code;
                }
                
                lmer_ &= lmer_mask_;
                
                // If we haven't filled first k-mer, don't return yet
                if ((str_pos_) >= (size_t)k_ && loaded_ch_ < l_) {
                    return &last_minimizer_;
                }
            }
            
            if (loaded_ch_ < l_) return nullptr;
            
            // Get canonical representation
            uint64_t canonical_lmer = canonical_representation(lmer_, l_);
            
            if (spaced_seed_mask_) {
                canonical_lmer &= spaced_seed_mask_;
            }
            
            uint64_t candidate_lmer = canonical_lmer ^ toggle_mask_;
            
            // Short-circuit for k == l
            if (k_ == l_) {
                last_minimizer_ = candidate_lmer ^ toggle_mask_;
                return &last_minimizer_;
            }
            
            // Sliding window minimum calculation
            while (!queue_.empty() && queue_.back().candidate > candidate_lmer) {
                queue_.pop_back();
            }
            
            MinimizerData data = {candidate_lmer, queue_pos_};
            
            if (queue_.empty() && queue_pos_ >= k_ - l_) {
                changed_minimizer = true;
            }
            
            queue_.push_back(data);
            
            // Expire l-mer not in current window
            if (!queue_.empty() && queue_.front().pos < queue_pos_ - k_ + l_) {
                queue_.pop_front();
                changed_minimizer = true;
            }
            
            // Change from no minimizer
            if (queue_pos_ == k_ - l_) {
                changed_minimizer = true;
            }
            
            queue_pos_++;
            
            if (str_pos_ >= (size_t)k_) {
                break;
            }
        }
        
        if (queue_.empty()) return nullptr;
        
        last_minimizer_ = queue_.front().candidate ^ toggle_mask_;
        return &last_minimizer_;
    }
    
    // Public accessors for parameters (needed for thread-local scanner creation)
    int k() const { return k_; }
    int l() const { return l_; }
    uint64_t spaced_seed_mask() const { return spaced_seed_mask_; }
    uint64_t toggle_mask() const { return toggle_mask_; }
    
private:
    struct MinimizerData {
        uint64_t candidate;
        int pos;
    };
    
    int k_;
    int l_;
    uint64_t spaced_seed_mask_;
    uint64_t toggle_mask_;
    uint64_t lmer_mask_;
    
    const std::string* str_;
    size_t str_pos_;
    size_t str_len_;
    int loaded_ch_;
    uint64_t lmer_;
    uint64_t last_minimizer_;
    
    std::deque<MinimizerData> queue_;
    int queue_pos_;
    
    uint8_t lookup_table_[256];
    
    void set_lookup(char c, uint8_t val) {
        lookup_table_[(unsigned char)c] = val;
    }
    
    // Kraken2's reverse complement implementation
    uint64_t reverse_complement(uint64_t kmer, int n) const {
        // Reverse bits (leaving bit pairs intact)
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCCULL) >> 2) | ((kmer & 0x3333333333333333ULL) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0ULL) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00ULL) >> 8) | ((kmer & 0x00FF00FF00FF00FFULL) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000ULL) >> 16) | ((kmer & 0x0000FFFF0000FFFFULL) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        
        // Complement and mask
        return ((~kmer) >> (64 - n * 2)) & ((1ULL << (n * 2)) - 1);
    }
    
    uint64_t canonical_representation(uint64_t kmer, int n) const {
        uint64_t revcom = reverse_complement(kmer, n);
        return kmer < revcom ? kmer : revcom;
    }
};

class FastaIndex {
public:
    struct IndexEntry {
        std::string name;
        size_t length;
        size_t offset;
        size_t line_bases;
        size_t line_width;
    };
    
    FastaIndex(const std::string& fasta_path) : fasta_path_(fasta_path) {
        std::string index_path = fasta_path + ".fai";
        
        if (!load_index(index_path)) {
            std::cout << "Creating FASTA index..." << std::endl;
            create_index();
            save_index(index_path);
        }
        
    }
    
    std::string get_sequence(const std::string& name) {
        auto it = index_.find(name);
        if (it == index_.end()) return "";
        
        const auto& entry = it->second;
        std::string sequence;
        sequence.reserve(entry.length);
        
        std::ifstream file(fasta_path_, std::ios::binary);
        if (!file) return "";
        
        file.seekg(entry.offset);
        
        size_t bases_read = 0;
        std::string line;
        while (bases_read < entry.length && std::getline(file, line)) {
            if (line.empty() || line[0] == '>') break;
            
            for (char c : line) {
                if (c != '\n' && c != '\r') {
                    sequence += c;
                    bases_read++;
                    if (bases_read >= entry.length) break;
                }
            }
        }
        
        return sequence;
    }
    
    size_t num_sequences() const { return index_.size(); }
    
private:
    std::string fasta_path_;
    std::unordered_map<std::string, IndexEntry> index_;
    
    void create_index() {
        std::ifstream file(fasta_path_);
        if (!file) throw std::runtime_error("Cannot open FASTA: " + fasta_path_);
        
        std::string line;
        std::string current_name;
        size_t current_offset = 0;
        size_t line_bases = 0;
        size_t line_width = 0;
        size_t seq_length = 0;
        size_t seq_start_offset = 0;
        
        while (std::getline(file, line)) {
            current_offset += line.length() + 1;
            
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (!current_name.empty()) {
                    IndexEntry entry;
                    entry.name = current_name;
                    entry.length = seq_length;
                    entry.offset = seq_start_offset;
                    entry.line_bases = line_bases;
                    entry.line_width = line_width;
                    index_[current_name] = entry;
                }
                
                current_name = line.substr(1);
                size_t space_pos = current_name.find(' ');
                if (space_pos != std::string::npos) {
                    current_name = current_name.substr(0, space_pos);
                }
                
                seq_length = 0;
                seq_start_offset = current_offset;
                line_bases = 0;
                line_width = 0;
            } else {
                size_t bases = 0;
                for (char c : line) {
                    if (c != '\n' && c != '\r') bases++;
                }
                
                seq_length += bases;
                if (line_bases == 0) {
                    line_bases = bases;
                    line_width = line.length() + 1;
                }
            }
        }
        
        if (!current_name.empty()) {
            IndexEntry entry;
            entry.name = current_name;
            entry.length = seq_length;
            entry.offset = seq_start_offset;
            entry.line_bases = line_bases;
            entry.line_width = line_width;
            index_[current_name] = entry;
        }
    }
    
    bool load_index(const std::string& index_path) {
        std::ifstream idx(index_path);
        if (!idx) return false;
        
        std::string line;
        while (std::getline(idx, line)) {
            std::istringstream iss(line);
            IndexEntry entry;
            iss >> entry.name >> entry.length >> entry.offset 
                >> entry.line_bases >> entry.line_width;
            index_[entry.name] = entry;
        }
        
        return !index_.empty();
    }
    
    void save_index(const std::string& index_path) {
        std::ofstream idx(index_path);
        for (const auto& [name, entry] : index_) {
            idx << entry.name << "\t" << entry.length << "\t" << entry.offset << "\t"
                << entry.line_bases << "\t" << entry.line_width << "\n";
        }
    }
};

class KrakenKmerExtractor {
public:
    struct ExtractedKmer {
        std::string sequence;
        uint32_t taxid;
    };
    
    KrakenKmerExtractor(int kmer_size, int minimizer_size, 
                        uint64_t spaced_seed_mask, uint64_t toggle_mask)
        : kmer_size_(kmer_size)
    {
        scanner_ = new MinimizerScanner(kmer_size, minimizer_size, 
                                        spaced_seed_mask, toggle_mask);
    }
    
    ~KrakenKmerExtractor() {
        delete scanner_;
    }
    
    std::vector<ExtractedKmer> extract_kmers(
        const std::string& kraken_file,
        FastaIndex& fasta_index,
        const TargetTaxIDManager& taxid_manager,
        int num_threads)
    {
        
        auto relevant_taxids = taxid_manager.get_all_relevant_taxids();
        
        size_t total_lines = count_lines(kraken_file);
        
        std::unordered_set<std::string> global_seen_kmers;
        std::vector<ExtractedKmer> all_kmers;
        std::mutex kmers_mutex;
        std::atomic<size_t> processed_count{0};
        
        const size_t BATCH_SIZE = 200000;
        
        std::ifstream kr_file(kraken_file);
        if (!kr_file) {
            throw std::runtime_error("Cannot open Kraken file: " + kraken_file);
        }
        
        std::vector<std::string> batch;
        batch.reserve(BATCH_SIZE);
        std::string line;
        size_t batch_num = 0;
        
        while (std::getline(kr_file, line)) {
            batch.push_back(std::move(line));
            
            if (batch.size() >= BATCH_SIZE) {
                process_batch(batch, fasta_index, taxid_manager, relevant_taxids, 
                             num_threads, global_seen_kmers, all_kmers, kmers_mutex, 
                             processed_count, batch_num++);
                batch.clear();
                
                std::cerr << "\rProcessed " << processed_count.load() << " / " 
                          << total_lines << " reads, extracted " 
                          << all_kmers.size() << " unique k-mers...";
            }
        }
        
        if (!batch.empty()) {
            process_batch(batch, fasta_index, taxid_manager, relevant_taxids,
                         num_threads, global_seen_kmers, all_kmers, kmers_mutex, 
                         processed_count, batch_num);
        }
        
        std::cerr << "\n";
        std::cout << "Extracted " << all_kmers.size() << " unique k-mers total" << std::endl;
        
        return all_kmers;
    }
    
private:
    int kmer_size_;
    MinimizerScanner* scanner_;
    
    void process_batch(
        const std::vector<std::string>& batch,
        FastaIndex& fasta_index,
        const TargetTaxIDManager& taxid_manager,
        const std::unordered_set<uint32_t>& relevant_taxids,
        int num_threads,
        std::unordered_set<std::string>& global_seen_kmers,
        std::vector<ExtractedKmer>& all_kmers,
        std::mutex& kmers_mutex,
        std::atomic<size_t>& processed_count,
        size_t /* batch_num */)
    {
        auto worker = [&](int tid) {
            std::vector<ExtractedKmer> local_kmers;
            std::unordered_set<std::string> local_seen;
            
            // Each thread needs its own scanner
            MinimizerScanner local_scanner(scanner_->k(), scanner_->l(), 
                                          scanner_->spaced_seed_mask(), 
                                          scanner_->toggle_mask());
            
            for (size_t i = tid; i < batch.size(); i += num_threads) {
                const auto& line = batch[i];
                
                std::istringstream iss(line);
                std::string classification, seq_id, taxid_str, seq_len_str, lca_mapping;
                
                if (!(iss >> classification >> seq_id >> taxid_str >> seq_len_str)) {
                    continue;
                }
                
                std::getline(iss, lca_mapping);
                if (!lca_mapping.empty() && lca_mapping[0] == '\t') {
                    lca_mapping = lca_mapping.substr(1);
                }
                
                std::string sequence = fasta_index.get_sequence(seq_id);
                if (sequence.empty()) continue;
                
                size_t seq_len = sequence.length();
                size_t kmer_read_index = 0;
                
                std::istringstream lca_stream(lca_mapping);
                std::string part;
                
                while (lca_stream >> part) {
                    size_t colon_pos = part.find(':');
                    if (colon_pos == std::string::npos) continue;
                    
                    std::string taxid_part = part.substr(0, colon_pos);
                    std::string count_part = part.substr(colon_pos + 1);
                    
                    if (taxid_part == "cov") continue;
                    
                    uint32_t taxid = 0;
                    int count = 0;
                    
                    try {
                        taxid = std::stoul(taxid_part);
                        count = std::stoi(count_part);
                    } catch (...) {
                        try {
                            count = std::stoi(count_part);
                        } catch (...) {
                            continue;
                        }
                        kmer_read_index += count;
                        continue;
                    }
                    
                    // Check if this taxid is relevant (either a target or descendant of rollup target)
                    if (relevant_taxids.find(taxid) != relevant_taxids.end()) {
                        // Get the target taxid to use (may be different due to rollup)
                        auto [should_include, target_taxid] = taxid_manager.get_target_for_taxid(taxid);
                        
                        if (should_include) {
                            for (int j = 0; j < count; ++j) {
                                size_t start_pos = kmer_read_index + j;
                                size_t end_pos = start_pos + kmer_size_;
                                
                                if (end_pos <= seq_len) {
                                    std::string kmer = sequence.substr(start_pos, kmer_size_);
                                    
                                    // Validate k-mer
                                    bool valid = true;
                                    for (char c : kmer) {
                                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                                            c != 'a' && c != 'c' && c != 'g' && c != 't') {
                                            valid = false;
                                            break;
                                        }
                                    }
                                    
                                    if (valid && local_seen.find(kmer) == local_seen.end()) {
                                        // Use the target_taxid (which may be rolled up)
                                        local_kmers.push_back({kmer, target_taxid});
                                        local_seen.insert(kmer);
                                    }
                                }
                            }
                        }
                    }
                    
                    kmer_read_index += count;
                }
                
                processed_count++;
            }
            
            std::lock_guard<std::mutex> lock(kmers_mutex);
            for (const auto& kmer : local_kmers) {
                if (global_seen_kmers.find(kmer.sequence) == global_seen_kmers.end()) {
                    all_kmers.push_back(kmer);
                    global_seen_kmers.insert(kmer.sequence);
                }
            }
        };
        
        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; ++i) {
            threads.emplace_back(worker, i);
        }
        for (auto& t : threads) t.join();
    }
    
    size_t count_lines(const std::string& filepath) {
        std::ifstream file(filepath);
        size_t count = 0;
        std::string line;
        while (std::getline(file, line)) count++;
        return count;
    }
};

class PackedTaxIDArray {
public:
    PackedTaxIDArray() = default;

    void build(size_t num_elements, int bits_per_element) {
        bits_per_element_ = bits_per_element;
        num_elements_ = num_elements;
        data_.resize((num_elements * bits_per_element + 63) / 64, 0);
    }

    void set(size_t index, uint64_t value) {
        if (index >= num_elements_) return;
        size_t bit_pos = index * bits_per_element_;
        size_t word_index = bit_pos / 64;
        size_t bit_offset = bit_pos % 64;
        
        uint64_t mask = (1ULL << bits_per_element_) - 1;
        value &= mask;

        data_[word_index] &= ~(mask << bit_offset);
        data_[word_index] |= (value << bit_offset);

        if (bit_offset + bits_per_element_ > 64 && word_index + 1 < data_.size()) {
            size_t remaining_bits = bit_offset + bits_per_element_ - 64;
            data_[word_index + 1] &= ~((1ULL << remaining_bits) - 1);
            data_[word_index + 1] |= (value >> (bits_per_element_ - remaining_bits));
        }
    }

    void save(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(&num_elements_), sizeof(num_elements_));
        out.write(reinterpret_cast<const char*>(&bits_per_element_), sizeof(bits_per_element_));
        size_t data_size = data_.size();
        out.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
        out.write(reinterpret_cast<const char*>(data_.data()), data_size * sizeof(uint64_t));
    }

private:
    std::vector<uint64_t> data_;
    size_t num_elements_ = 0;
    int bits_per_element_ = 0;
};

class KmerDatabaseBuilder {
public:
    KmerDatabaseBuilder(int k, int l, uint64_t spaced_seed_mask, uint64_t toggle_mask,
                        int tax_id_bits, int fingerprint_bits);
    ~KmerDatabaseBuilder();
    
    void build_from_kmers(const std::vector<KrakenKmerExtractor::ExtractedKmer>& kmers, 
                          int num_threads);
    void save_to_disk(const std::string& db_prefix) const;
    
private:
    using MPHF_Hasher = boomphf::SingleHashFunctor<kmer_t>;
    using MPHF_Type = boomphf::mphf<kmer_t, MPHF_Hasher>;

    int k_;
    int l_;
    uint64_t spaced_seed_mask_;
    uint64_t toggle_mask_;
    int tax_id_bits_;
    int fingerprint_bits_;
    
    MPHF_Type* mphf_;
    std::vector<uint16_t> fingerprint_array_;
    PackedTaxIDArray tax_id_array_;
    size_t num_unique_minimizers_;
    
    std::vector<uint32_t> index_to_taxid_;
    std::unordered_map<uint32_t, uint32_t> taxid_to_index_;
    
    MinimizerScanner* scanner_;
    
    uint16_t get_fingerprint(kmer_t kmer) const;
};

KmerDatabaseBuilder::KmerDatabaseBuilder(int k, int l, uint64_t spaced_seed_mask, 
                                         uint64_t toggle_mask, int tax_id_bits, 
                                         int fingerprint_bits)
    : k_(k), l_(l), spaced_seed_mask_(spaced_seed_mask), toggle_mask_(toggle_mask),
      tax_id_bits_(tax_id_bits), fingerprint_bits_(fingerprint_bits), 
      mphf_(nullptr), num_unique_minimizers_(0)
{
    if (l <= 0 || l > 31) throw std::invalid_argument("l must be 1-31");
    if (k < l) throw std::invalid_argument("k must be >= l");
    if (fingerprint_bits <= 0 || fingerprint_bits > 16) 
        throw std::invalid_argument("Fingerprint bits must be 1-16");
    
    scanner_ = new MinimizerScanner(k, l, spaced_seed_mask, toggle_mask);
}

KmerDatabaseBuilder::~KmerDatabaseBuilder() {
    delete mphf_;
    delete scanner_;
}

void KmerDatabaseBuilder::build_from_kmers(
    const std::vector<KrakenKmerExtractor::ExtractedKmer>& kmers, 
    int num_threads) 
{
    std::cout << "Input: " << kmers.size() << " k-mers" << std::endl;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // Phase 1: TaxID mapping
    std::cout << "\n[Phase 1/4] Creating TaxID mapping..." << std::endl;
    std::set<uint32_t> unique_taxids;
    for (const auto& kmer : kmers) {
        unique_taxids.insert(kmer.taxid);
    }
        
    if (unique_taxids.size() > (1ULL << tax_id_bits_)) {
        throw std::runtime_error("Too many unique taxids for " + std::to_string(tax_id_bits_) + " bits");
    }
    
    uint32_t index = 0;
    for (uint32_t taxid : unique_taxids) {
        index_to_taxid_.push_back(taxid);
        taxid_to_index_[taxid] = index;
        std::cout << "  TaxID " << taxid << " -> index " << index << std::endl;
        index++;
    }
        
    std::unordered_map<kmer_t, uint32_t> minimizer_map;
    std::mutex map_mutex;
    
    auto worker = [&](int tid) {
        std::unordered_map<kmer_t, uint32_t> local_map;
        MinimizerScanner local_scanner(k_, l_, spaced_seed_mask_, toggle_mask_);
        
        for (size_t i = tid; i < kmers.size(); i += num_threads) {
            local_scanner.load_sequence(kmers[i].sequence);
            
            uint64_t* minimizer_ptr = local_scanner.next_minimizer();
            if (minimizer_ptr) {
                uint32_t tax_id_index = taxid_to_index_[kmers[i].taxid];
                local_map.insert({*minimizer_ptr, tax_id_index});
            }
        }
        
        std::lock_guard<std::mutex> lock(map_mutex);
        for (const auto& [k, v] : local_map) {
            minimizer_map.insert({k, v});
        }
    };
    
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(worker, i);
    }
    for (auto& t : threads) t.join();
    
    num_unique_minimizers_ = minimizer_map.size();
    std::cout << "  Extracted " << num_unique_minimizers_ << " unique minimizers" << std::endl;
    
    if (num_unique_minimizers_ == 0) {
        throw std::runtime_error("No minimizers found!");
    }
    
    
    std::vector<kmer_t> keys;
    keys.reserve(num_unique_minimizers_);
    for (const auto& [k, _] : minimizer_map) {
        keys.push_back(k);
    }
    
    mphf_ = new MPHF_Type(num_unique_minimizers_, keys, num_threads, 2.0);
        
    fingerprint_array_.resize(num_unique_minimizers_);
    tax_id_array_.build(num_unique_minimizers_, tax_id_bits_);
    
    std::vector<std::pair<kmer_t, uint32_t>> map_items(minimizer_map.begin(), minimizer_map.end());
    
    auto populate_worker = [&](int tid) {
        for (size_t i = tid; i < map_items.size(); i += num_threads) {
            const auto& [kmer, tax_id_index] = map_items[i];
            uint64_t idx = mphf_->lookup(kmer);
            if (idx < num_unique_minimizers_) {
                fingerprint_array_[idx] = get_fingerprint(kmer);
                tax_id_array_.set(idx, tax_id_index);
            }
        }
    };
    
    std::vector<std::thread> populate_threads;
    for (int i = 0; i < num_threads; ++i) {
        populate_threads.emplace_back(populate_worker, i);
    }
    for (auto& t : populate_threads) t.join();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(t_end - t_start).count();
    
    std::cout << "Total time: " << duration << " seconds" << std::endl;
    std::cout << "Unique minimizers: " << num_unique_minimizers_ << std::endl;
}

void KmerDatabaseBuilder::save_to_disk(const std::string& db_prefix) const {
    if (!mphf_) throw std::runtime_error("Cannot save empty database");
        
    std::ofstream meta(db_prefix + ".meta", std::ios::binary);
    meta.write(reinterpret_cast<const char*>(&k_), sizeof(k_));
    meta.write(reinterpret_cast<const char*>(&l_), sizeof(l_));
    meta.write(reinterpret_cast<const char*>(&spaced_seed_mask_), sizeof(spaced_seed_mask_));
    meta.write(reinterpret_cast<const char*>(&toggle_mask_), sizeof(toggle_mask_));
    meta.write(reinterpret_cast<const char*>(&tax_id_bits_), sizeof(tax_id_bits_));
    meta.write(reinterpret_cast<const char*>(&fingerprint_bits_), sizeof(fingerprint_bits_));
    meta.write(reinterpret_cast<const char*>(&num_unique_minimizers_), sizeof(num_unique_minimizers_));
    meta.close();
    
    std::ofstream mphf_file(db_prefix + ".mphf", std::ios::binary);
    mphf_->save(mphf_file);
    mphf_file.close();
    
    std::ofstream fp(db_prefix + ".fp", std::ios::binary);
    fp.write(reinterpret_cast<const char*>(fingerprint_array_.data()), 
             fingerprint_array_.size() * sizeof(uint16_t));
    fp.close();
    
    std::ofstream taxid(db_prefix + ".taxid", std::ios::binary);
    tax_id_array_.save(taxid);
    taxid.close();
    
    std::ofstream taxmap(db_prefix + ".taxmap", std::ios::binary);
    size_t map_size = index_to_taxid_.size();
    taxmap.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    taxmap.write(reinterpret_cast<const char*>(index_to_taxid_.data()), 
                 map_size * sizeof(uint32_t));
    taxmap.close();
    
    std::ofstream taxmap_txt(db_prefix + ".taxmap.txt");
    taxmap_txt << "Index\tActual_TaxID\n";
    for (size_t i = 0; i < index_to_taxid_.size(); ++i) {
        taxmap_txt << i << "\t" << index_to_taxid_[i] << "\n";
    }
    taxmap_txt.close();
    
}

uint16_t KmerDatabaseBuilder::get_fingerprint(kmer_t kmer) const {
    kmer ^= kmer >> 33;
    kmer *= 0xff51afd7ed558ccdULL;
    kmer ^= kmer >> 33;
    return static_cast<uint16_t>(kmer & ((1ULL << fingerprint_bits_) - 1));
}

uint64_t create_spaced_seed_mask_from_pattern(const std::string& bit_pattern) {
    int l = bit_pattern.length();
    uint64_t mask = 0;
    
    for (int pos = 0; pos < l; ++pos) {
        char bit = bit_pattern[l - 1 - pos];
        if (bit == '1') {
            mask |= (3ULL << (2 * pos));
        }
    }
    
    return mask;
}

std::unordered_set<uint32_t> load_target_taxids(const std::string& target_file) {
    std::unordered_set<uint32_t> targets;
    
    std::ifstream file(target_file);
    if (!file) {
        throw std::runtime_error("Cannot open target file: " + target_file);
    }
    
    std::string line;
    std::getline(file, line); // Skip header
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        uint32_t taxid;
        if (iss >> taxid) {
            if (taxid > 1) {
                targets.insert(taxid);
            }
        }
    }
    
    std::cout << "Loaded " << targets.size() << " target taxids" << std::endl;
    return targets;
}

void print_usage(const char* prog_name) {
    std::cerr << "K-mer Database Builder\n\n"
              << "Usage:\n"
              << "  " << prog_name << " <kraken_output> <fasta_file> <targets_file> <nodes_dmp> <db_prefix> [threads]\n\n"
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string kraken_file = argv[1];
    std::string fasta_file = argv[2];
    std::string targets_file = argv[3];
    std::string nodes_file = argv[4];
    std::string db_prefix = argv[5];
    int threads = (argc >= 7) ? std::atoi(argv[6]) : std::thread::hardware_concurrency();
    
    const int K = 35;
    const int L = 31;
    const int TAXID_BITS = 8;
    const int FP_BITS = 16;
    
    std::string spaced_pattern = "1111111111111111111110101010101";
    uint64_t spaced_seed_mask = create_spaced_seed_mask_from_pattern(spaced_pattern);
    uint64_t toggle_mask = 0;
    
    
    std::cout << std::endl;
    
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        TaxonomyTree taxonomy;
        taxonomy.load_from_nodes_dmp(nodes_file);
        
        auto targets = load_target_taxids(targets_file);
        
        TargetTaxIDManager taxid_manager(targets, taxonomy);
        
        FastaIndex fasta_index(fasta_file);
        
        KrakenKmerExtractor extractor(K, L, spaced_seed_mask, toggle_mask);
        auto kmers = extractor.extract_kmers(kraken_file, fasta_index, taxid_manager, threads);
        
        KmerDatabaseBuilder builder(K, L, spaced_seed_mask, toggle_mask, TAXID_BITS, FP_BITS);
        builder.build_from_kmers(kmers, threads);
        builder.save_to_disk(db_prefix);
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration<double>(total_end - total_start).count();
        
        
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}