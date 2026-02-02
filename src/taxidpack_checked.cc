class PackedTaxIDArray {
public:
//I think I have to (for every subtype species node/taxid) also sample all of its children. 
//Make sure all those taxids are independent.

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

    uint64_t get(size_t index) const {
        if (index >= num_elements_) return 0;
        size_t bit_pos = index * bits_per_element_;
        size_t word_index = bit_pos / 64;
        size_t bit_offset = bit_pos % 64;

        uint64_t mask = (1ULL << bits_per_element_) - 1;
        uint64_t value = (data_[word_index] >> bit_offset);

        if (bit_offset + bits_per_element_ > 64 && word_index + 1 < data_.size()) {
            size_t remaining_bits = bit_offset + bits_per_element_ - 64;
            value |= (data_[word_index + 1] & ((1ULL << remaining_bits) - 1)) 
                     << (bits_per_element_ - remaining_bits);
        }
        return value & mask;
    }

    void save(std::ofstream& out) const {

        // I think writing size_t directly on disk might cause some vulnerability since diff machiines might have 32bit
        out.write(reinterpret_cast<const char*>(&num_elements_), sizeof(num_elements_));
        out.write(reinterpret_cast<const char*>(&bits_per_element_), sizeof(bits_per_element_));
        size_t data_size = data_.size();
        out.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
        out.write(reinterpret_cast<const char*>(data_.data()), data_size * sizeof(uint64_t));
    }

    void load(std::ifstream& in) {
        in.read(reinterpret_cast<char*>(&num_elements_), sizeof(num_elements_));
        in.read(reinterpret_cast<char*>(&bits_per_element_), sizeof(bits_per_element_));
        size_t data_size;
        in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
        data_.resize(data_size);
        in.read(reinterpret_cast<char*>(data_.data()), data_size * sizeof(uint64_t));
    }

private:
    int target_count_;
    std::vector<uint64_t> data_;
    size_t num_elements_ = 0;
    int bits_per_element_ = 6;
};