/*
 * Code to build the SA and LCP arrays of the dictionary of a prefix-free parse.
 * 
 * This code is adapted from https://github.com/maxrossi91/pfp-thresholds/blob/master/include/pfp/dictionary.hpp
 */

#ifndef DICTIONARY_HPP
#define DICTIONARY_HPP

#include "pfp_parse.hpp"

extern "C" {
    #include "gsacak.h"
}

class dictionary{
private:
    std::vector<uint32_t> fchar;
    std::vector<bool> word_is_sequence_start; // Flag for each word in the dictionary

    
    
    
public:
    std::vector<uint8_t> d;
    std::vector<uint32_t> occ; // Algo limit, a word cannot occurs more than 2^32 times
    std::vector<uint_t> saD;
    std::vector<int_t> lcpD; // Int because of gsacak interface
    sdsl::bit_vector b_d; // Starting position of each phrase in D    
    sdsl::bit_vector b_s;
    sdsl::bit_vector::rank_1_type rank_b_d;
    sdsl::bit_vector::select_1_type select_b_d;
    sdsl::bit_vector::rank_1_type rank_b_s;

    // function that returns true if the suffix is valid false otherwise
    inline bool is_valid_suffix(size_t pos, size_t len) {
        // Check for terminators in the suffix
        for (size_t i = 0; i < len; i++) {
            if (pos + i < d.size() && d[pos + i] == '$') {
                return false; // Found a terminator in the suffix
            }
        }
        return true; // No terminator found
    }

    bool isWordSequenceStart(size_t word_index) const {
        if (word_index < word_is_sequence_start.size()) {
            return word_is_sequence_start[word_index];
        }
        return false;
    }

    // default constructor for load.
    dictionary() {}

    dictionary(std::string filename, size_t w) 
    {
        std::string tmp_filename = filename + std::string(".edict");
        read_file(tmp_filename.c_str(), d);
        std::cout << "Dictionary loaded: " << d.size() << " bytes" << std::endl;
        
        tmp_filename = filename + std::string(".eocc");
        read_file(tmp_filename.c_str(), occ);
        std::cout << "Occurrences loaded: " << occ.size() << " elements" << std::endl;
        
        tmp_filename = filename + std::string(".fchar");
        read_file(tmp_filename.c_str(), fchar);
        std::cout << "First chars loaded: " << fchar.size() << " elements" << std::endl;
        
        std::cout << "Building bitvector..." << std::endl;

        // Bitvector for sentence beginnings (b_d)
        b_d = sdsl::bit_vector(d.size(), 0);
        b_d[0] = 1;  // First sentence always starts at 0

        for(size_t i = 1; i < d.size() - 1; i++) {
            if(d[i] == EndOfWord) {
                b_d[i + 1] = 1;  // Next sentence starts after EndOfWord
            }
        }

        // Bitvector for sequence beginnings (b_s)
        b_s = sdsl::bit_vector(d.size(), 0);
        
        rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
        select_b_d = sdsl::bit_vector::select_1_type(&b_d);
        
       // Configure b_s based on fchar data
       for(size_t i = 0; i < fchar.size(); i += 2) {
            size_t pos = select_b_d(fchar[i]) + fchar[i+1] +1;
            if (pos < d.size()) {
                b_s[pos] = 1;
                
                std::string word_at_pos = "";
                for (size_t j = 0; j < 20 && pos + j < d.size(); j++) {
                    if (d[pos + j] == EndOfWord || d[pos + j] == EndOfDict) break;
                    word_at_pos += static_cast<char>(d[pos + j]);
                }
            }
        }
        
        rank_b_s = sdsl::bit_vector::rank_1_type(&b_s);
        
        std::string seqstart_filename = filename + ".seqstart";
        FILE* seqstart_file = fopen(seqstart_filename.c_str(), "rb");
        if (seqstart_file) {
            size_t num_starts;
            if (fread(&num_starts, sizeof(num_starts), 1, seqstart_file) == 1) {
                std::vector<uint64_t> seq_start_hashes(num_starts);
                for (size_t i = 0; i < num_starts; i++) {
                    if (fread(&seq_start_hashes[i], sizeof(uint64_t), 1, seqstart_file) != 1) {
                        break;
                    }
                }
            }
            fclose(seqstart_file);
        } 
        
        word_is_sequence_start.resize(occ.size(), false);
        
        
        // Build SA and LCP for the original dictionary
        build();
        
       
    }

    void build(){
        std::cout << "Building dictionary SA and LCP arrays..." << std::endl;
        
        saD.resize(d.size());
        lcpD.resize(d.size());
        
        verbose("Computing SA and LCP of the dictionary");
        _elapsed_time(
            std::cout << "Starting gsacak for dictionary size: " << d.size() << std::endl;
            gsacak(&d[0], &saD[0], &lcpD[0], nullptr, d.size());
            std::cout << "gsacak completed" << std::endl;
        );
        
    }
};

#endif /* DICTIONARY_HPP */