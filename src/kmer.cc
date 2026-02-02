class KmerDatabase {
public:
    KmerDatabase(int k, int l, uint64_t spaced_seed_mask, uint64_t toggle_mask, int fingerprint_bits);
    ~KmerDatabase() { delete mphf_; }
    
    void build_from_fasta(const std::string& fasta_path, int num_threads);
    void save_to_disk(const std::string& db_prefix) const;
    bool load_from_disk(const std::string& db_prefix);
    std::map<uint32_t, int> query_read(const std::string& read_sequence) const;
    uint32_t get_actual_taxid(uint32_t index) const;
    
private:
    using MPHF_Hasher = boomphf::SingleHashFunctor<kmer_t>;
    using MPHF_Type = boomphf::mphf<kmer_t, MPHF_Hasher>;

    int k_;
    int l_;
    uint64_t spaced_seed_mask_;
    uint64_t toggle_mask_;
    int tax_id_bits_ = 
    int fingerprint_bits_;
    uint64_t lmer_mask_;
    
    MPHF_Type* mphf_ = nullptr;
    std::vector<uint16_t> fingerprint_array_;
    PackedTaxIDArray tax_id_array_;
    size_t num_unique_minimizers_ = 0;
    
    // TaxID mapping: index -> actual taxid
    std::vector<uint32_t> index_to_taxid_;
    std::unordered_map<uint32_t, uint32_t> taxid_to_index_;
    
    uint64_t reverse_complement(uint64_t kmer, int n) const;
    uint64_t canonical_representation(uint64_t kmer, int n) const;
    uint16_t get_fingerprint(kmer_t kmer) const;
    uint64_t* get_minimizer_from_sequence(const std::string& sequence, 
                                          uint64_t& minimizer_out) const;
    std::vector<kmer_t> get_minimizers_from_read(const std::string& read) const;
};

KmerDatabase::KmerDatabase(int k, int l, uint64_t spaced_seed_mask, uint64_t toggle_mask,
                           int tax_id_bits, int fingerprint_bits)
    : k_(k), l_(l), spaced_seed_mask_(spaced_seed_mask), toggle_mask_(toggle_mask),
      tax_id_bits_(tax_id_bits), fingerprint_bits_(fingerprint_bits)
{
    if (l <= 0 || l > 31) throw std::invalid_argument("l must be 1-31");
    if (k < l) throw std::invalid_argument("k must be >= l");
    if (fingerprint_bits <= 0 || fingerprint_bits > 16) 
        throw std::invalid_argument("Fingerprint bits must be 1-16");
    
    lmer_mask_ = (1ULL << (l_ * 2)) - 1;
    toggle_mask_ &= lmer_mask_;
}

std::vector<kmer_t> KmerDatabase::get_minimizers_from_read(const std::string& read) const {
    std::vector<kmer_t> minimizers;
    size_t read_len = read.length();
    if (read_len < k_) return minimizers;

    std::deque<std::pair<uint64_t, int>> window_queue; // <hash, l-mer_start_pos>
    uint64_t lmer = 0;
    int valid_bases = 0;
    
    size_t i = 0; // Current base position (end of base)
    
    
    #if defined(__AVX2__)
        // Use 32-byte AVX2 registers
        const __m256i low_nibble_mask = _mm256_set1_epi8(0x0F);
        const __m256i lut = _mm256_setr_epi8(
            4, 0, 4, 1, 3, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 16 bytes
            4, 0, 4, 1, 3, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4  // 16 bytes
        );
        alignas(32) uint8_t codes[32];
        const int SIMD_WIDTH = 32;
    #elif defined(__SSE3__) // pshufb is SSE3
        // Use 16-byte SSE3 registers
        const __m128i low_nibble_mask = _mm_set1_epi8(0x0F);
        const __m128i lut = _mm_setr_epi8(
            4, 0, 4, 1, 3, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4
        );
        alignas(16) uint8_t codes[16];
        const int SIMD_WIDTH = 16;
    #elif defined(__ARM_NEON) || defined(__ARM_NEON__)
        // Use 16-byte ARM NEON registers
        const uint8x16_t low_nibble_mask = vdupq_n_u8(0x0F);
        const uint8_t lut_array[16] = {
            4, 0, 4, 1, 3, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4
        };
        const uint8x16_t lut = vld1q_u8(lut_array);
        alignas(16) uint8_t codes[16];
        const int SIMD_WIDTH = 16;
    #else
        // Scalar fallback if no SIMD is available/enabled
        const int SIMD_WIDTH = 1;
        uint8_t codes[1]; // Dummy array
    #endif

    // --- Main loop ---
    // This structure processes in chunks of SIMD_WIDTH, then handles
    // the remainder. If SIMD_WIDTH=1 (scalar), it just runs as a
    // normal for loop.
    size_t main_len = read_len - (read_len % SIMD_WIDTH);
    for (i = 0; i < main_len; i += SIMD_WIDTH) {
        #if defined(__AVX2__)
            __m256i chars = _mm256_loadu_si256((__m256i*)(&read[i]));
            __m256i nibbles = _mm256_and_si256(chars, low_nibble_mask);
            __m256i simd_codes = _mm256_shuffle_epi8(lut, nibbles);
            _mm256_store_si256((__m256i*)codes, simd_codes);
        #elif defined(__SSE3__)
            __m128i chars = _mm_loadu_si128((__m128i*)(&read[i]));
            __m128i nibbles = _mm_and_si128(chars, low_nibble_mask);
            __m128i simd_codes = _mm_shuffle_epi8(lut, nibbles);
            _mm_store_si128((__m128i*)codes, simd_codes);
        #elif defined(__ARM_NEON) || defined(__ARM_NEON__)
            uint8x16_t chars = vld1q_u8(reinterpret_cast<const uint8_t*>(&read[i]));
            uint8x16_t nibbles = vandq_u8(chars, low_nibble_mask);
            uint8x16_t simd_codes = vqtbl1q_u8(lut, nibbles);
            vst1q_u8(codes, simd_codes);
        #endif

        // Serial part: process the converted codes
        // This loop is the $O(N)$ sliding window logic
        for (int j = 0; j < SIMD_WIDTH; ++j) {
            uint8_t code;
            int current_pos = i + j;

            #if defined(__AVX2__) || defined(__SSE3__) || defined(__ARM_NEON) || defined(__ARM_NEON__)
                code = codes[j];
            #else
                // Scalar fallback path
                char c = read[current_pos];
                switch(c) {
                    case 'A': case 'a': code = 0; break;
                    case 'C': case 'c': code = 1; break;
                    case 'G': case 'g': code = 2; break;
                    case 'T': case 't': code = 3; break;
                    default: code = 4;
                }
            #endif

            if (code > 3) { // Handle 'N' or other invalid
                lmer = 0;
                valid_bases = 0;
                window_queue.clear();
                continue;
            }

            lmer = ((lmer << 2) | code) & lmer_mask_;
            valid_bases++;

            if (valid_bases < l_) continue; // Wait for the first full l-mer

            // We have a valid l-mer ending at `current_pos`.
            // It *starts* at `current_pos - l + 1`.
            int lmer_start_pos = current_pos - l_ + 1;
            
            uint64_t canonical = canonical_representation(lmer, l_);
            if (spaced_seed_mask_) {
                canonical &= spaced_seed_mask_;
            }
            uint64_t candidate = canonical ^ toggle_mask_;
            
            // Add new l-mer to deque
            while (!window_queue.empty() && window_queue.back().first > candidate) {
                window_queue.pop_back();
            }
            window_queue.push_back({candidate, lmer_start_pos});
            
            // A k-mer *ending* at `current_pos` starts at `current_pos - k + 1`.
            // The l-mers in this k-mer window start from `[current_pos - k + 1]`
            // up to `[current_pos - l + 1]`.
            int window_start_pos = current_pos - k_ + 1;
            
            // Remove l-mers that are no longer in this k-mer's window
            while (!window_queue.empty() && window_queue.front().second < window_start_pos) {
                window_queue.pop_front();
            }
            
            // Have we processed a full k-mer?
            // The first k-mer *ends* at position `k - 1`.
            if (current_pos >= k_ - 1) {
                // The front of the queue is the minimizer for the k-mer
                // *ending* at `current_pos`.
                minimizers.push_back(window_queue.front().first ^ toggle_mask_);
            }
        }
    }
    
    // --- Scalar epilogue for remaining bases ---
    // This loop handles the few bases left over if read_len % SIMD_WIDTH != 0
    // If we used the scalar path, i == read_len, so this loop is skipped.
    for (; i < read_len; ++i) {
        uint8_t code;
        char c = read[i];
        switch(c) {
            case 'A': case 'a': code = 0; break;
            case 'C': case 'c': code = 1; break;
            case 'G': case 'g': code = 2; break;
            case 'T': case 't': code = 3; break;
            default:
                lmer = 0;
                valid_bases = 0;
                window_queue.clear();
                continue;
        }
        
        lmer = ((lmer << 2) | code) & lmer_mask_;
        valid_bases++;

        if (valid_bases < l_) continue;

        int lmer_start_pos = i - l_ + 1;
        
        uint64_t canonical = canonical_representation(lmer, l_);
        if (spaced_seed_mask_) {
            canonical &= spaced_seed_mask_;
        }
        uint64_t candidate = canonical ^ toggle_mask_;
        
        while (!window_queue.empty() && window_queue.back().first > candidate) {
            window_queue.pop_back();
        }
        window_queue.push_back({candidate, lmer_start_pos});
        
        int window_start_pos = i - k_ + 1;
        while (!window_queue.empty() && window_queue.front().second < window_start_pos) {
            window_queue.pop_front();
        }
        
        if (i >= k_ - 1) {
            minimizers.push_back(window_queue.front().first ^ toggle_mask_);
        }
    }

    return minimizers;
}
// === NEW CODE END ===

void KmerDatabase::build_from_fasta(const std::string& fasta_path, int num_threads) {
    std::cout << "Reading " << k_ << "-mers from " << fasta_path << "..." << std::endl;
    
    std::ifstream infile(fasta_path);
    if (!infile) throw std::runtime_error("Cannot open file: " + fasta_path);
    
    std::vector<std::pair<std::string, uint32_t>> sequences;
    std::set<uint32_t> unique_taxids;  // Track unique taxids
    std::string line, header;
    
    while (std::getline(infile, header)) {
        if (header.empty() || header[0] != '>') continue;
        if (!std::getline(infile, line)) break;
        
        // Parse header: >NC_075450.1|kmer_start_0|taxid_10509
        size_t last_pipe = header.rfind('|');
        if (last_pipe == std::string::npos) continue;
        
        std::string taxid_part = header.substr(last_pipe + 1);
        if (taxid_part.find("taxid_") != 0) continue;
        
        uint32_t tax_id = std::stoul(taxid_part.substr(6));
        unique_taxids.insert(tax_id);
        
        if ((int)line.length() != k_) continue;
        
        bool valid = true;
        for (char c : line) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                valid = false;
                break;
            }
        }
        
        if (valid) sequences.push_back({line, tax_id});
    }
    
    std::cout << "Number of k-mers:" << sequences.size() << std::endl;
    std::cout << "Number of unique taxids: " << unique_taxids.size() << std::endl;
    
    // Create TaxID mapping: assign sequential indices
    uint32_t index = 0;
    for (uint32_t taxid : unique_taxids) {
        index_to_taxid_.push_back(taxid);
        taxid_to_index_[taxid] = index;
        std::cout << "  TaxID " << taxid << " -> index " << index << std::endl;
        index++;
    }
    
    if (unique_taxids.size() > (1ULL << 64)) {
        throw std::runtime_error("There is no way. You are living in a world with aliens! " + std::to_string(tax_id_bits_) + " bits");
    }
    
    std::unordered_map<kmer_t, uint32_t> minimizer_map;
    std::mutex map_mutex;
    
    auto worker = [&](int tid) {
        std::unordered_map<kmer_t, uint32_t> local_map;
        for (size_t i = tid; i < sequences.size(); i += num_threads) {
            uint64_t minimizer;
            // NOTE: The build path still uses the O(K) function, which is correct
            // because it operates on one k-mer at a time.
            if (get_minimizer_from_sequence(sequences[i].first, minimizer)) { 
                // Convert actual taxid to index
                uint32_t tax_id_index = taxid_to_index_[sequences[i].second];
                local_map.insert({minimizer, tax_id_index});
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
    std::cout << "Extracted " << num_unique_minimizers_ << " unique minimizers" << std::endl;
    
    if (num_unique_minimizers_ == 0) {
        throw std::runtime_error("No minimizers found!");
    }
    
    std::vector<kmer_t> keys;
    keys.reserve(num_unique_minimizers_);
    for (const auto& [k, _] : minimizer_map) keys.push_back(k);
    
    std::cout << "Building MPHF..." << std::endl;
    mphf_ = new MPHF_Type(num_unique_minimizers_, keys, num_threads, 2.0);
    
    std::cout << "Populating arrays..." << std::endl;
    fingerprint_array_.resize(num_unique_minimizers_);
    tax_id_array_.build(num_unique_minimizers_, tax_id_bits_);
    
    for (const auto& [kmer, tax_id_index] : minimizer_map) {
        uint64_t idx = mphf_->lookup(kmer);
        if (idx < num_unique_minimizers_) {
            fingerprint_array_[idx] = get_fingerprint(kmer);
            tax_id_array_.set(idx, tax_id_index);
        }
    }
    
    std::cout << "Database built successfully!" << std::endl;
}

void KmerDatabase::save_to_disk(const std::string& db_prefix) const {
    if (!mphf_) throw std::runtime_error("Cannot save empty database");
    
    std::cout << "Saving database to disk with prefix: " << db_prefix << std::endl;
    
    // Save metadata
    std::ofstream meta(db_prefix + ".meta", std::ios::binary);
    meta.write(reinterpret_cast<const char*>(&k_), sizeof(k_));
    meta.write(reinterpret_cast<const char*>(&l_), sizeof(l_));
    meta.write(reinterpret_cast<const char*>(&spaced_seed_mask_), sizeof(spaced_seed_mask_));
    meta.write(reinterpret_cast<const char*>(&toggle_mask_), sizeof(toggle_mask_));
    meta.write(reinterpret_cast<const char*>(&tax_id_bits_), sizeof(tax_id_bits_));
    meta.write(reinterpret_cast<const char*>(&fingerprint_bits_), sizeof(fingerprint_bits_));
    meta.write(reinterpret_cast<const char*>(&num_unique_minimizers_), sizeof(num_unique_minimizers_));
    meta.close();
    
    // Save MPHF
    std::ofstream mphf_file(db_prefix + ".mphf", std::ios::binary);
    mphf_->save(mphf_file);
    mphf_file.close();
    
    // Save fingerprints
    std::ofstream fp(db_prefix + ".fp", std::ios::binary);
    fp.write(reinterpret_cast<const char*>(fingerprint_array_.data()), 
             fingerprint_array_.size() * sizeof(uint16_t));
    fp.close();
    
    // Save taxids
    std::ofstream taxid(db_prefix + ".taxid", std::ios::binary);
    tax_id_array_.save(taxid);
    taxid.close();
    
    // Save taxid mapping
    std::ofstream taxmap(db_prefix + ".taxmap", std::ios::binary);
    size_t map_size = index_to_taxid_.size();
    taxmap.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    taxmap.write(reinterpret_cast<const char*>(index_to_taxid_.data()), 
                 map_size * sizeof(uint32_t));
    taxmap.close();
    
    // Save human-readable taxid mapping
    std::ofstream taxmap_txt(db_prefix + ".taxmap.txt");
    taxmap_txt << "Index\tActual_TaxID\n";
    for (size_t i = 0; i < index_to_taxid_.size(); ++i) {
        taxmap_txt << i << "\t" << index_to_taxid_[i] << "\n";
    }
    taxmap_txt.close();
    
    std::cout << "Database saved successfully!" << std::endl;
    std::cout << "  Files: " << db_prefix << ".{meta,mphf,fp,taxid,taxmap,taxmap.txt}" << std::endl;
}

bool KmerDatabase::load_from_disk(const std::string& db_prefix) {
    std::cout << "Loading database from disk: " << db_prefix << std::endl;
    
    try {
        // Load metadata
        std::ifstream meta(db_prefix + ".meta", std::ios::binary);
        if (!meta) return false;
        meta.read(reinterpret_cast<char*>(&k_), sizeof(k_));
        meta.read(reinterpret_cast<char*>(&l_), sizeof(l_));
        meta.read(reinterpret_cast<char*>(&spaced_seed_mask_), sizeof(spaced_seed_mask_));
        meta.read(reinterpret_cast<char*>(&toggle_mask_), sizeof(toggle_mask_));
        meta.read(reinterpret_cast<char*>(&tax_id_bits_), sizeof(tax_id_bits_));
        meta.read(reinterpret_cast<char*>(&fingerprint_bits_), sizeof(fingerprint_bits_));
        meta.read(reinterpret_cast<char*>(&num_unique_minimizers_), sizeof(num_unique_minimizers_));
        meta.close();
        
        lmer_mask_ = (1ULL << (l_ * 2)) - 1;
        toggle_mask_ &= lmer_mask_;
        
        // Load MPHF
        std::ifstream mphf_file(db_prefix + ".mphf", std::ios::binary);
        if (!mphf_file) return false;
        delete mphf_;
        mphf_ = new MPHF_Type();
        mphf_->load(mphf_file);
        mphf_file.close();
        
        // Load fingerprints
        std::ifstream fp(db_prefix + ".fp", std::ios::binary);
        if (!fp) return false;
        fingerprint_array_.resize(num_unique_minimizers_);
        fp.read(reinterpret_cast<char*>(fingerprint_array_.data()), 
                fingerprint_array_.size() * sizeof(uint16_t));
        fp.close();
        
        // Load taxids
        std::ifstream taxid(db_prefix + ".taxid", std::ios::binary);
        if (!taxid) return false;
        tax_id_array_.load(taxid);
        taxid.close();
        
        // Load taxid mapping
        std::ifstream taxmap(db_prefix + ".taxmap", std::ios::binary);
        if (!taxmap) return false;
        size_t map_size;
        taxmap.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        index_to_taxid_.resize(map_size);
        taxmap.read(reinterpret_cast<char*>(index_to_taxid_.data()), 
                    map_size * sizeof(uint32_t));
        taxmap.close();
        
        // Rebuild reverse map
        taxid_to_index_.clear();
        for (size_t i = 0; i < index_to_taxid_.size(); ++i) {
            taxid_to_index_[index_to_taxid_[i]] = i;
        }
        
        std::cout << "Database loaded successfully!" << std::endl;
        std::cout << "  K=" << k_ << ", L=" << l_ << ", Unique minimizers=" 
                  << num_unique_minimizers_ << std::endl;
        std::cout << "  TaxID mapping: " << index_to_taxid_.size() << " unique taxa" << std::endl;
        return true;
    } catch (...) {
        return false;
    }
}

uint32_t KmerDatabase::get_actual_taxid(uint32_t index) const {
    if (index < index_to_taxid_.size()) {
        return index_to_taxid_[index];
    }
    return 0;
}

std::map<uint32_t, int> KmerDatabase::query_read(const std::string& read_sequence) const {
    std::map<uint32_t, int> counts;
    if (!mphf_ || (int)read_sequence.length() < k_) return counts;
    
    std::vector<kmer_t> query_minimizers = get_minimizers_from_read(read_sequence);
    
    for (kmer_t minimizer : query_minimizers) {
        uint64_t idx = mphf_->lookup(minimizer);
        
        if (idx < num_unique_minimizers_) {
            if (fingerprint_array_[idx] == get_fingerprint(minimizer)) {
                uint32_t tax_id = tax_id_array_.get(idx);
                counts[tax_id]++;
            }
        }
    }
    
    return counts;
}

uint16_t KmerDatabase::get_fingerprint(kmer_t kmer) const {
    kmer ^= kmer >> 33;
    kmer *= 0xff51afd7ed558ccdULL;
    kmer ^= kmer >> 33;
    return static_cast<uint16_t>(kmer & ((1ULL << fingerprint_bits_) - 1));
}