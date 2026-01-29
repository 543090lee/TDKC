#include <zlib.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread)

class FastqReader {
public:
    FastqReader(const std::string& filename) : gz_file_(nullptr), ks_(nullptr) {
        //gzopen should technically handle both .gz and plain
        gz_file_ = gzopen(filename.c_str(), "r");
        if (!gz_file_) {
            throw std::runtime_error("Can't open your file man");
        }
        ks_ = kseq_init(gz_file_);
    }

    ~FastqReader() {
        if (ks_) kseq_destroy(ks_);
        if (gz_file_) gzclose(gz_file_);
    }

    bool next_read(std::string& header, std::string& sequence, std::string& quality) {
        //greater than zero should be a success
        if (kseq_read(ks_) < 0) {
            return false;
        }

        //name.s should hold the header ID, @ to end, and comment after
        if (ks_->comment.l) {
            header.assign(ks_->name.s);
            header += " ";
            header += ks_->comment.s;
        } else {
            header.assign(ks_->name.s);
        }

        sequence.assign(ks_->seq.s, ks_->seq.l);

        if (ks_->qual.l) {
            quality.assign(ks_->qual.s, ks_->qual.l);
        } else {
            quality.clear();
        }

        return true;
    }

private:
    gzFile gz_file_;
    //ks parser from kseq
    kseq_t* ks_; 
};