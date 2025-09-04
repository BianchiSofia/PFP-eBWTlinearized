#ifndef PFP_HPP
#define PFP_HPP

#define BWTBYTES 5

#include <algorithm>
#include <tuple>
#include <queue>
#include <cctype>
#include <string>
#include <limits>
#include <set>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <cstdlib>


// Enum per i tipi di posizioni
enum PositionType {
    TYPE_L = 0,
    TYPE_S = 1,
    TYPE_LMS = 2
};

// Forward declarations
class pfp_parse;
class dictionary;

class pfp {
private:
    typedef struct
    {
        size_t i = 0;                        // Index in the suffix array
        size_t phrase = 0;
        size_t suffix_length = 0;
        size_t sn = 0;                       // Position in the dictionary
        uint8_t bwt_char = 0;                // Changed from std::string to uint8_t
        size_t st_pos = 0;
        bool is_starting = 0;
        PositionType type = TYPE_S;          // Position type (S, L, LMS)
    } phrase_suffix_t;
    
    // Vectors to store the types of the positions
    std::vector<PositionType> position_types;
    std::vector<bool> is_lms_position;

    std::vector<std::string> debug_words;
    FILE* debug_suffixes_file = nullptr;
    
    size_t first_occ = 0; // First occurrence of same suffix phrases in BWT_P
    size_t last_occ = 0;     // Last occurrence of same suffix phrases in BWT_P
    size_t text_len = 0;
    size_t ins_sofar = 0; // Number of characters inserted in the eBWT
    size_t runs = 0; // Number of runs of the eBWT
    
    bool rle;
    
    // for rle eBWT
    FILE *ebwt_file_len;
    FILE *ebwt_file_heads;
    // for plain eBWT
    FILE *ebwt_file;
    FILE *I_file;

    void calculatePositionTypes() {
        position_types.resize(dict.d.size(), TYPE_S);
        is_lms_position.resize(dict.d.size(), false);
        
        for (int i = dict.d.size() - 2; i >= 0; i--) {
            if (dict.d[i] < dict.d[i + 1]) {
                position_types[i] = TYPE_S;
            } else if (dict.d[i] > dict.d[i + 1]) {
                position_types[i] = TYPE_L;
            } else {
                position_types[i] = position_types[i + 1];
            }
        }
        
        for (size_t i = 1; i < dict.d.size(); i++) {
            if (position_types[i] == TYPE_S && position_types[i - 1] == TYPE_L) {
                is_lms_position[i] = true;
                position_types[i] = TYPE_LMS;
            }
        }
    }

   
public:
    pfp_parse& pars;
    dictionary& dict;
    std::string filename;
    size_t w;
    size_t min_s; // Value of the minimum lcp_T in the current run of BWT_T
    size_t pos_s; // Position of the minimum lcp_T in the current run of BWT_T

    uint8_t head; 
    size_t start;
    size_t length = 0; // Length of the current run of BWT_T
    
    // Check if a suffix is valid 
    inline bool is_valid(phrase_suffix_t& s) {
        for (size_t i = 0; i < std::min(s.suffix_length, size_t(20)) && s.sn + i < dict.d.size(); i++) {
            if (dict.d[s.sn + i] == 1 || dict.d[s.sn + i] == 0) { 
                break;
            }
        }
        if (s.suffix_length == 0) {
            return false;
        }
        
        bool is_phrase_start = dict.b_d[s.sn] == 1; 
        if (is_phrase_start) {
            return false;
        }
        
        bool contains_dollar = false;
        for (size_t i = 0; i < s.suffix_length && s.sn + i < dict.d.size(); i++) {
            uint8_t char_val = dict.d[s.sn + i];
            if (char_val == 36) { 
                contains_dollar = true;
                break;
            } 
        }
        
        if (contains_dollar) {
            return true;
        }

        if (s.is_starting) {
            if (s.suffix_length < w) {
                return false;
            }
        } else {
            if (s.suffix_length < w) {
                return false;
            }
        }

        return true;
    }
    
    inline void update_ebwt(uint8_t next_char, size_t length_, bool is_st, size_t ind, size_t st_pos)
    {
        if (head != next_char)
        {
            print_ebwt();

            head = next_char;
            length = 0; // Debug only
        }
        
        ins_sofar += length_;
        length += length_;
        
        if(is_st && st_pos > 0){
            if(pars.b_st[ind]==1){
                if(pars.offset[(pars.rank_st(ind+1)-1)] == st_pos){
                    start = ins_sofar-1;
                    if(fwrite(&start,sizeof(start),1,I_file)!=1) 
                        throw std::runtime_error("I file write error");
                }
            }
        }
    }

    // Increment to the next suffix 
    inline bool inc(phrase_suffix_t& s) {
        s.i++;
        if (s.i >= dict.saD.size())
            return false;
        
        s.sn = dict.saD[s.i];
        s.is_starting = dict.b_s[s.sn];
        bool is_phrase_start = dict.b_d[s.sn] == 1;
        s.phrase = dict.rank_b_d(s.sn+1);
            
        if(s.is_starting) { 
            s.st_pos = s.sn - dict.select_b_d(s.phrase); 
        } else { 
            s.st_pos = 0; 
        }
        
        s.suffix_length = dict.select_b_d(dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        
        assert(s.phrase > 0 && s.phrase <= pars.ilP.size());
        
        // Assign position type
        if (s.sn < position_types.size()) {
            s.type = position_types[s.sn];
        }
        
        if(is_valid(s)) {
            if (s.sn > 0) {
                uint8_t prev_char = dict.d[s.sn - 1];
                
                s.bwt_char = prev_char;
            } 
        } else {
            s.bwt_char = 0;
        }
        
        return true;
    }
    
    // Print eBWT run 
    inline void print_ebwt() {
        if(length > 0)
        {
            if(rle)
            {
                // Write the head
                if (fputc(head, ebwt_file_heads) == EOF)
                    throw std::runtime_error("BWT write error 1");
                
                // Write the length
                if (fwrite(&length, BWTBYTES, 1, ebwt_file_len) != 1)
                    throw std::runtime_error("BWT write error 2");

            }else{
                // if not rle write plain ebwt
                for(size_t i = 0; i < length; ++i)
                {
                    if (fputc(head, ebwt_file) == EOF)
                        throw std::runtime_error("BWT write error 1");
                }
            }
            // one run added
            ++runs;
        }
    }

    pfp(pfp_parse &p_, dictionary &d_, std::string filename_, size_t w_, bool rle_): 
        pars(p_),
        dict(d_),
        filename(filename_),
        w(w_),
        rle(rle_),
        pos_s(0),
        head(0)  
    {
        assert(dict.d[dict.saD[0]] == 0); 
        calculatePositionTypes();

        std::string outfile = filename + std::string(".I");
        if((I_file = fopen(outfile.c_str(), "w")) == nullptr)
            throw std::runtime_error("open() file " + outfile + " failed"); // Usa throw invece di error()
        
        if(rle)
        {    
            outfile = filename + std::string(".ebwt.heads");
            if ((ebwt_file_heads = fopen(outfile.c_str(), "w")) == nullptr)
                throw std::runtime_error("open() file " + outfile + " failed");

            outfile = filename + std::string(".ebwt.len");
            if ((ebwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                throw std::runtime_error("open() file " + outfile + " failed");
        }else{

            outfile = filename + std::string(".ebwt");
            if ((ebwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                throw std::runtime_error("open() file " + outfile + " failed");
        }

        phrase_suffix_t curr;
        
        inc(curr);
        
        size_t valid_suffixes = 0;
        
        while (curr.i < dict.saD.size())
        {
            if(is_valid(curr))
            {
                valid_suffixes++;
                size_t start_pos = curr.sn;
                size_t max_len = std::min(curr.suffix_length, (size_t)50);
                size_t printed = 0;

                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);  // Store the list of all phrase ids with the same suffix.

                bool same_chars = true;
                bool st_chars = curr.is_starting;

                phrase_suffix_t next = curr;

                size_t collected = 1;
            
                while (inc(next) && (dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((dict.b_d[next.sn] == 0 && next.suffix_length >= w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_chars = (same_chars && same_suffix.back().bwt_char == next.bwt_char);
                        st_chars = (st_chars || next.is_starting);
                        same_suffix.push_back(next);                       
                    }                  
                }

                for (size_t i = 0; i < same_suffix.size(); i++) {
                    const auto& s = same_suffix[i];
                }

                // Simple case
                if (same_chars && !st_chars){
                    for (auto curr : same_suffix)
                    {
                        update_ebwt(curr.bwt_char, dict.occ[curr.phrase-1], 0, 0, 0);
                    }
                }   
                // Hard case
                else
                {       
                    //suffix not starting with a character occurring at the beginning of 
                    // a input sequence
                    if(!st_chars){
                        for (size_t i = 0; i < same_suffix.size(); i++) {
                            auto s = same_suffix[i];
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                        }

                        typedef std::pair<uint_s *, std::pair<uint_s *, uint8_t>> pq_t;
                        // using lambda to compare elements.
                        auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                            return *lhs.first > *rhs.first;
                        };
                        std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                        
                        for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            pq.push(std::make_pair(&pars.ilP[begin], std::make_pair(&pars.ilP[end], s.bwt_char)));
                        }
                            
                        while(!pq.empty()){
                            auto curr_occ = pq.top();
                            pq.pop();
                            update_ebwt(curr_occ.second.second, 1, 0, 0, 0);
                            curr_occ.first++;
                            if (curr_occ.first != curr_occ.second.first ){ 
                                pq.push(curr_occ);
                            }
                        }
                    } else {
                        //check and store the positions at which the start of string characters
                        //are inserted in the ebwt
                        typedef std::pair<uint_s *, std::tuple<uint_s *, uint8_t, size_t, size_t>> tq_t;
                        

                        auto cmp2 = [](const tq_t &lhs, const tq_t &rhs) {
                            return *lhs.first > *rhs.first;
                        };
                        std::priority_queue<tq_t, std::vector<tq_t>, decltype(cmp2)> tq(cmp2);
                        
                        
                        for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            tq.push(std::make_pair(&pars.ilP[begin], std::make_tuple(&pars.ilP[end], s.bwt_char, begin, s.st_pos)));
                        }
                        
                        while(!tq.empty()){
                            auto curr_occ = tq.top();
                            tq.pop();
                            
                            update_ebwt(std::get<1>(curr_occ.second), 1, 1, std::get<2>(curr_occ.second), std::get<3>(curr_occ.second));
                            // Update pq
                            curr_occ.first++;
                            std::get<2>(curr_occ.second)++;
                            if (curr_occ.first != std::get<0>(curr_occ.second)){
                                tq.push(curr_occ);
                            }
                        
                        }
                    }
                }
                curr = next;
            }
            else
            {
                inc(curr);
            }      
        }
        // print last run
        print_ebwt();
        // close output files
        fclose(I_file);
        // close rle lengths file if rle was used
        if(rle){ fclose(ebwt_file_len); fclose(ebwt_file_heads);}
        else { fclose(ebwt_file); }   

        std::string info = "Number of runs of the eBWT: " + std::to_string(runs) + "\n"; 
        info += "Length of the eBWT: " + std::to_string(ins_sofar) + "\n"; 
        if(rle)
        {
            info += "Number of bytes per eBWT run character (.heads file): 1\n";
            info += "Number of bytes per eBWT run length (.len file): " + std::to_string(BWTBYTES) + "\n";
        } 
        else
            info += "Number of bytes per eBWT character (.ebwt file): 1\n";

        std::ofstream out(filename + std::string(".info"));
        out << info;
        out.close();
        
    }      
};

#endif /* PFP_HPP */