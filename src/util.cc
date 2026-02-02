
uint64_t fit_bit_pattern_to_64_ssm(const std::string& bit_pattern) {
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


uint64_t reverse_complement(uint64_t kmer, int n) const {
    kmer = ((kmer & 0xCCCCCCCCCCCCCCCCULL) >> 2) | ((kmer & 0x3333333333333333ULL) << 2);
    kmer = ((kmer & 0xF0F0F0F0F0F0F0F0ULL) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4);
    kmer = ((kmer & 0xFF00FF00FF00FF00ULL) >> 8) | ((kmer & 0x00FF00FF00FF00FFULL) << 8);
    kmer = ((kmer & 0xFFFF0000FFFF0000ULL) >> 16) | ((kmer & 0x0000FFFF0000FFFFULL) << 16);
    kmer = (kmer >> 32) | (kmer << 32);
    return ((~kmer) >> (64 - n * 2)) & ((1ULL << (n * 2)) - 1);
}

uint64_t canonical_representation(uint64_t kmer, int n) const {
    uint64_t revcom = reverse_complement(kmer, n);
    return kmer < revcom ? kmer : revcom;
}

void print_usage_querying(const char* prog_name) {
    std::cerr << "Usage:\n"
              << "  Query reads:\n"
              << "    " << prog_name << " query <db_prefix> <your fastq[.gz]> [threads]\n"
}

void print_usage_building(const char* prog_name) {
    std::cerr << "Usage:\n"
              << "  Build db:\n"
              << "    " << prog_name << " build <input.fasta> <self-kraken.output> <taxa-list.tsv> <db_prefix> [threads]\n"
}