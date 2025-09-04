/*
 * PFP parse implementation to compute the non-circular Prefix-free parse of a collection of sequences.
 * 
 * This code is adapted from https://github.com/alshai/Big-BWT/blob/master/newscan.cpp
 *
 */
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <map>
#include <set>
#include <assert.h>
#include <errno.h>
#include <zlib.h>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
#include "xerrors.h"
}
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;
typedef pair <uint32_t,uint32_t> p;

// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   size_t w = 10;            // sliding window size and its default
   size_t p = 100;           // modulus for establishing stopping w-tuples
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level
};

struct word_stats {
  string str;  // parse phrase
  occ_int_t occ;  // no. of phrases
  word_int_t rank=0; // rank of the phrase
};

// bitvector to keep track of the positions of the '$' terminators
struct BitvectorTerminators {
  vector<uint64_t> positions; 
  vector<int> sequence_numbers;

  // Modify the function to match the call site
  void add_terminator_position(uint64_t pos, int seq_num = -1) {
      positions.push_back(pos);
      sequence_numbers.push_back(seq_num);
  }
  
  void save_to_file(const string& filename) {
      FILE *f = fopen((filename + ".term").c_str(), "wb");
      if (!f) die("Cannot open terminator positions file for writing");
      
      uint64_t size = positions.size();
      if (fwrite(&size, sizeof(size), 1, f) != 1) 
          die("Error writing terminator count");
      
      for (uint64_t pos : positions) {
          if (fwrite(&pos, sizeof(pos), 1, f) != 1)
              die("Error writing terminator position");
      }
      
      fclose(f);
  }
};

// Global bitvector 
BitvectorTerminators terminatorPositions;
std::map<std::string, int> terminator_sequence_order; 
int global_sequence_counter = 0; 


void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        #ifndef NOTHREADS
        << "\t-t M\tnumber of helper threads, def. none " << endl
        #endif
        << "\t-h  \tshow help and exit" << endl;
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
    int c;
    extern char *optarg;
    extern int optind;

    puts("==== Command line:");
    for(int i=0;i<argc;i++)
      printf(" %s",argv[i]);
    puts("");
  
    string sarg;
    while ((c = getopt( argc, argv, "p:w:ht:v") ) != -1) {
       switch(c) {
         case 'w':
         sarg.assign( optarg );
         arg.w = stoi( sarg ); break;
         case 'p':
         sarg.assign( optarg );
         arg.p = stoi( sarg ); break;
         case 't':
         sarg.assign( optarg );
         arg.th = stoi( sarg ); break;
         case 'v':
            arg.verbose++; break;
         case 'h':
            print_help(argv, arg); exit(1);
         case '?':
         cout << "Unknown option. Use -h for help." << endl;
         exit(1);
       }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
      arg.inputFileName.assign( argv[optind] );
    }
    else {
       cout << "Invalid number of arguments" << endl;
       print_help(argv,arg);
    }
    // check algorithm parameters
    if(arg.w <4) {
      cout << "Windows size must be at least 4\n";
      exit(1);
    }
    if(arg.p<10) {
      cout << "Modulus must be at leas 10\n";
      exit(1);
    }
    #ifdef NOTHREADS
    if(arg.th!=0) {
      cout << "The NT version cannot use threads\n";
      exit(1);
    }
    #else
    if(arg.th<0) {
      cout << "Number of threads cannot be negative\n";
      exit(1);
    }
    #endif
}

struct KR_window {
  int wsize;
  int current;
  int *window;
  int asize;
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime

  KR_window(int w): wsize(w) {
    asize = 256;
    asize_pot = 1;
    for(int i=1;i<wsize;i++)
      asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for(int i=0;i<wsize;i++) window[i]=0;
    // init hash value and related values
    hash=tot_char=current=0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    current++;
    current = min(wsize,current);
    // complex expression to avoid negative numbers
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
    hash = (asize*hash + c) % prime;      //  add char i
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash;
  }
  // debug only
  string get_window() {
    string w = "";
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }

  ~KR_window() {
    delete[] window;
  }
};

static void save_update_word(string& w, unsigned int minsize, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, bool last_word);

#ifndef NOTHREADS
#include "circpfp.hpp"
#endif

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      assert(c>=0 && c< 256);
      hash = (256*hash + c) % prime;    //  add char k
    }
    return hash;
}

static void save_update_word(string& w, unsigned int minsize, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, bool last_word)
{
  if(w.size() <= minsize) {
    return;
  }
  
  // Get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  
  if(fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1) 
    die("parse write error");
  
  // Update frequency table for current hash
  if (freq.find(hash) == freq.end()) {
    freq[hash].occ = 1; // new hash
    freq[hash].str = w;
  } else {
    freq[hash].occ += 1; // known hash - + occ
    
    if (freq[hash].occ <= 0) {
      cerr << "Emergency exit! Maximum # of occurrence of dictionary word exceeded\n";
      exit(1);
    }
    if (freq[hash].str != w) {
      cerr << "Emergency exit! Hash collision for strings:\n";
      cerr << freq[hash].str << "\n  vs\n" <<  w << endl;
      exit(1);
    }
  }
  
  if (last_word) {
    return;
  }
  
  // keep only the overlapping part of the window
  w.erase(0, w.size() - minsize);
}




// Non-circular prefix free parse of fname, w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t parse_fasta(Args& arg, map<uint64_t,word_stats>& wordFreq)
{
    string fnam = arg.inputFileName;
  
    // open the 1st pass parsing file
    FILE *parse_file = open_aux_file(arg.inputFileName.c_str(),EXTPARS0,"wb");
    // open the words offset file
    FILE *offset_file = open_aux_file(arg.inputFileName.c_str(),EXTOFF0,"wb");
  
    // main loop on the chars of the input file
    uint8_t c;
    uint64_t total_char = 0;
    uint64_t sequence_pos = 0;
    uint64_t seq_number = 0;
    
    vector<size_t> sequence_start_hashes;
    
    gzFile fp;
    kseq_t *seq;
    long int l;
    fp = gzopen(fnam.c_str(), "r");
    seq = kseq_init(fp);
    int seqn=0;
    
    while ((l = kseq_read(seq)) >= 0) {
        seqn++;
        seq_number++;
        
        string current_sequence(seq->seq.s, seq->seq.l);
        
        KR_window krw(arg.w); 
        string current_word = ""; 
        bool found_first_trigger = false;
        uint64_t last_trigger_pos = 0;
        
        uint64_t sequence_absolute_start = 0;
        if (fwrite(&sequence_absolute_start, sizeof(sequence_absolute_start), 1, offset_file) != 1) 
            die("offset write error");
        
        // Process each character in the sequence
        for (size_t i = 0; i < seq->seq.l; i++) {
            c = std::toupper(seq->seq.s[i]);
            if (c <= Dollar) {
                cerr << "Invalid char found in input file: no additional chars will be read\n";
                break;
            }
            
            current_word.append(1, c);
            uint64_t hash = krw.addchar(c);
            
            // Check if this is a trigger position (hash % p == 0 and window is full)
            if (hash % arg.p == 0 && krw.current == arg.w) {
                
                if (!found_first_trigger) {
                    // First trigger in this sequence
                    found_first_trigger = true;
                    
                    current_word.insert(0, "$"); // Add terminator at the start
  
                    // Add the first w-tuple to the dictionary
                    uint64_t first_hash = kr_hash(current_word);
                    if (fwrite(&first_hash, sizeof(first_hash), 1, parse_file) != 1) 
                        die("parse write error");

              
                    sequence_start_hashes.push_back(first_hash);
                    
                    // Update frequency table
                    if (wordFreq.find(first_hash) == wordFreq.end()) {
                        wordFreq[first_hash].occ = 1;
                        wordFreq[first_hash].str = current_word;
              
                    } else {
                        wordFreq[first_hash].occ += 1;
                        if (wordFreq[first_hash].str != current_word) {
                            cerr << "Hash collision detected!\n";
                            exit(1);
                        }
                    }
                    
                    // Keep only the last w characters for the next phrase
                    last_trigger_pos = i;
                    current_word = current_word.substr(current_word.size() - arg.w);
                                
                  } else {
                    // Not the first trigger, save the word from last trigger to here
              
                    save_update_word(current_word, arg.w, wordFreq, parse_file, false);
                    
                    // Keep only the last w characters for the next phrase
                    last_trigger_pos = i;
                    current_word = current_word.substr(current_word.size() - arg.w);

                }
            }
        }
        
        if (!found_first_trigger) {
            string terminator = "$";
            current_word.append(terminator);
           
            terminator_sequence_order[current_word] = global_sequence_counter++;
            current_word.insert(0, terminator); // Add terminator at the start

            // Add to dictionary
            uint64_t hash = kr_hash(current_word);
            if (fwrite(&hash, sizeof(hash), 1, parse_file) != 1) 
                die("parse write error");
            
            sequence_start_hashes.push_back(hash);
            
            // Update frequency table
            if (wordFreq.find(hash) == wordFreq.end()) {
                wordFreq[hash].occ = 1;
                wordFreq[hash].str = current_word;
            } else {
                wordFreq[hash].occ += 1;
                if (wordFreq[hash].str != current_word) {
                    cerr << "Hash collision detected!\n";
                    exit(1);
                }
            }
        } else {
            string terminator = "$";
            current_word.append(terminator);
            
            terminator_sequence_order[current_word] = global_sequence_counter++;
            
            // Save the final word to the dictionary 
            save_update_word(current_word, arg.w, wordFreq, parse_file, true);
        }
        
        // Save the position of the terminator
        terminatorPositions.add_terminator_position(sequence_pos + krw.tot_char, seq_number);
        
        total_char += krw.tot_char;
        sequence_pos += krw.tot_char;

    }
    
    kseq_destroy(seq);
    gzclose(fp);

    // close input and output files
    if(fclose(parse_file)!=0) die("Error closing parse file");
    if(fclose(offset_file)!=0) die("Error closing offset file");
    
    // Save the bitvector of terminator positions
    terminatorPositions.save_to_file(arg.inputFileName);
    
    std::string seqstart_filename = arg.inputFileName + ".seqstart";
    FILE* seqstart_file = fopen(seqstart_filename.c_str(), "wb");
    if (!seqstart_file) die("Cannot open sequence start file for writing");
    
    size_t num_starts = sequence_start_hashes.size();
    if (fwrite(&num_starts, sizeof(num_starts), 1, seqstart_file) != 1)
        die("Error writing number of sequence starts");
    
    for (uint64_t hash : sequence_start_hashes) {
        if (fwrite(&hash, sizeof(hash), 1, seqstart_file) != 1)
            die("Error writing sequence start hash");
    }
    
    fclose(seqstart_file);
    
    return total_char;
}

bool pstringCompare(const string *a, const string *b);

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict)
{  
  vector<pair<string, word_stats*>> dict_entries;
  for (auto& entry : wfreq) {
    dict_entries.push_back({entry.second.str, &entry.second});
  }
  
  sort(dict_entries.begin(), dict_entries.end(), 
       [](const auto& a, const auto& b) {
         return pstringCompare(&a.first, &b.first);
       });
  
  FILE *fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICT, "wb");
  FILE *focc = open_aux_file(arg.inputFileName.c_str(), EXTOCC, "wb");

  word_int_t wrank = 1; // current word rank (1 based)
  
  for (const auto& entry : dict_entries) {
    const string& word = entry.first;
    word_stats* wf = entry.second;
    
    size_t len = word.size();
    assert(len > (size_t)arg.w);
    
    size_t s = fwrite(word.data(), 1, len, fdict);
    if(s != len) die("Error writing to DICT file");
    if(fputc(EndOfWord, fdict) == EOF) die("Error writing EndOfWord to DICT file");
    
    assert(wf->occ > 0);
    s = fwrite(&wf->occ, sizeof(wf->occ), 1, focc);
    if(s != 1) die("Error writing to OCC file");
    
    assert(wf->rank == 0);  
    wf->rank = wrank++;
  }
  
  if(fputc(EndOfDict, fdict) == EOF) die("Error writing EndOfDict to DICT file");
  if(fclose(focc) != 0) die("Error closing OCC file");
  if(fclose(fdict) != 0) die("Error closing DICT file");
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
    const string& sa = *a;
    const string& sb = *b;
    
    for (size_t i = 0; i < min(sa.size(), sb.size()); i++) {
        if (sa[i] == '$' && sb[i] != '$') return true;
        if (sa[i] != '$' && sb[i] == '$') return false;
        
        if (sa[i] == '$' && sb[i] == '$') {
            auto it_a = terminator_sequence_order.find(sa);
            auto it_b = terminator_sequence_order.find(sb);
            
            if (it_a != terminator_sequence_order.end() && 
                it_b != terminator_sequence_order.end()) {
                return it_a->second < it_b->second;
            }
            
            return sa < sb;  // lex
        }
        
        if (sa[i] != sb[i]) return sa[i] < sb[i];
    }
    
    return sa.size() < sb.size();
}

void remapParse(Args &arg, map<uint64_t,word_stats> &wfreq, int th)
{  
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, th);
  mFile *moff = mopen_aux_file(arg.inputFileName.c_str(), EXTOFF0, th);
  FILE *newp   = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");
  FILE *newoff  = open_aux_file(arg.inputFileName.c_str(), EXTOFF, "wb");
  FILE *strt   = open_aux_file(arg.inputFileName.c_str(), EXTSTART, "wb");
  FILE *fchar  = open_aux_file(arg.inputFileName.c_str(), EXTFCHAR, "wb");
  
  // recompute occ as an extra check
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based
  uint64_t hash, phash, fc;
  uint64_t start = 0, len = 0;
  set<p> startChr;
  
  vector<uint64_t> offsets;
  while(true) {
    uint64_t offset;
    size_t s = mfread(&offset, sizeof(offset), 1, moff);
    if(s == 0) break;
    if(s != 1) die("Error reading offset file");
    offsets.push_back(offset);
  }
  
  if(mfclose(moff) != 0) die("Error closing offset file");
  moff = mopen_aux_file(arg.inputFileName.c_str(), EXTOFF0, th);
  
  size_t offset_idx = 0;
  bool expecting_sequence_start = true; 
  size_t word_counter = 0;
  
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    
    word_counter++;
    len++;
    
    auto it = wfreq.find(hash);
    if (it == wfreq.end()) {
      die("Hash not found in updated dictionary");
    }
    
    word_int_t rank = it->second.rank;
    string current_word = it->second.str;
  
    bool contains_dollar = (current_word.find('$') != string::npos);
    bool starts_with_dollar = (!current_word.empty() && current_word[0] == '$');
    bool ends_with_dollar = (!current_word.empty() && current_word.back() == '$');
      
    occ[rank]++;
    phash = hash;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to new parse file");
    
    if (expecting_sequence_start) {
        s = fwrite(&start,sizeof(start),1,strt);
        if(s!=1) die("Error writing to start file");
        p st = p(rank, 0);
        if(startChr.find(st) == startChr.end()) {
            startChr.insert(st);
        }
        
        expecting_sequence_start = false;
    }
    
    if (ends_with_dollar) {
        if(offset_idx < offsets.size()) {
            fc = offsets[offset_idx];
            size_t last_dollar = current_word.rfind('$');
            uint32_t off = (last_dollar != string::npos) ? uint32_t(last_dollar) : uint32_t(fc);
            
            s = fwrite(&off,sizeof(off),1,newoff);
            if(s!=1) die("Error writing to new offset file");
            
            offset_idx++;
        } else {
            uint32_t off = 0;
            s = fwrite(&off,sizeof(off),1,newoff);
            if(s!=1) die("Error writing to new offset file");
        }
        
        start += len;
        len = 0;
        
        expecting_sequence_start = true;
    }
  }

  for (auto& x: startChr) {
    if(fwrite(&x,sizeof(x),1,fchar)!=1) die("error writing to first char file");
  }
  
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(fclose(fchar)!=0) die("Error closing first char positions file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  if(mfclose(moff)!=0) die("Error closing offset file");
  if(fclose(strt)!=0) die("Error closing starting positions file");
  if(fclose(newoff)!=0) die("Error closing new offsets file");
}

int main(int argc, char** argv) {
    
    // translate command line parameters
    Args arg;
    parseArgs(argc, argv, arg);

    // measure elapsed wall clock time
    time_t start_main = time(NULL);
    time_t start_wc = start_main;
    // init sorted map counting the number of occurrences of parse phrases
    map <uint64_t,word_stats> wordFreq;
    uint64_t totChar; int nt = 0; // tot characters seen
    
    // ------------ parse input fasta file
    try{
        if(arg.th<=1){totChar = parse_fasta(arg,wordFreq);}
        else
        {
            #ifdef NOTHREADS
            cerr << "Sorry, this is the no-threads executable and you requested " << arg.th << " threads\n";
            exit(1);
            #else
            Res res = parallel_parse_fasta(arg, wordFreq);
            totChar = res.tot_char; nt = res.us_th;
            for (uint64_t pos : res.terminator_positions) {
              terminatorPositions.add_terminator_position(pos);
            }
            #endif
        }
    }
    catch(const std::bad_alloc&) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
    }
    
    uint64_t totDWord = wordFreq.size();
    cout << "Total input symbols: " << totChar << endl;
    cout << "Found " << totDWord << " distinct words" <<endl;
    cout << "Number of sequences: " << terminatorPositions.positions.size() << endl;
    cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    // check # distinct words
    if(totDWord>MAX_DISTINCT_WORDS) {
      cerr << "Emergency exit! The number of distinc words (" << totDWord << ")\n";
      cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
      exit(1);
    }
    
    // -------------- second pass
    start_wc = time(NULL);
    // create array of dictionary words
    vector<const string *> dictArray;
    dictArray.reserve(totDWord);
    // fill array
    uint64_t sumLen = 0;
    uint64_t totWord = 0;
    for (auto& x: wordFreq) {
      sumLen += x.second.str.size();
      totWord += x.second.occ;
      dictArray.push_back(&x.second.str);
    }
    assert(dictArray.size()==totDWord);
    cout << "Sum of lenghts of dictionary words: " << sumLen << endl;
    cout << "Total number of words: " << totWord << endl;
    
    vector<string> sorted_phrases;

    for (auto& x: wordFreq) {
      size_t dollar_pos = x.second.str.find('$');
      if (dollar_pos != string::npos) {
          string word = x.second.str;
          
          if (dollar_pos == word.size() - 1) {
              sorted_phrases.push_back({word}); 
          } else {
              sorted_phrases.push_back({word});
          }
      }
    }
    sort(sorted_phrases.begin(), sorted_phrases.end(), 
        [](const string& a, const string& b) { return pstringCompare(&a, &b); });

    sort(dictArray.begin(), dictArray.end(), pstringCompare);
    
    for (size_t i = 0; i < dictArray.size(); i++) {
        std::string word = *dictArray[i];
        size_t dollar_pos = word.find('$');
    }

    // write plain dictionary and occ file, also compute rank for each hash
    cout << "Writing plain dictionary and occ file\n";
    writeDictOcc(arg, wordFreq, dictArray);
    dictArray.clear(); // reclaim memory
    cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    

    std::map<int, uint64_t> original_terminator_positions; // seq_num -> posizione
    for (size_t i = 0; i < terminatorPositions.positions.size(); i++) {
        int seq_num = i + 1;
        uint64_t pos = terminatorPositions.positions[i];
        original_terminator_positions[seq_num] = pos;
    }

    std::string termpos_filename = arg.inputFileName + ".termpos";
    FILE* termpos_file = fopen(termpos_filename.c_str(), "wb");
    if (!termpos_file) die("Cannot open terminator positions file for writing");

    size_t num_terms = original_terminator_positions.size();
    if (fwrite(&num_terms, sizeof(num_terms), 1, termpos_file) != 1)
        die("Error writing number of terminators");

    for (const auto& [seq_num, pos] : original_terminator_positions) {
        if (fwrite(&seq_num, sizeof(seq_num), 1, termpos_file) != 1)
            die("Error writing sequence number");
        if (fwrite(&pos, sizeof(pos), 1, termpos_file) != 1)
            die("Error writing terminator position");
    }
    fclose(termpos_file);
    
    // remap parse file
    start_wc = time(NULL);
    cout << "Generating remapped parse file\n";
    remapParse(arg, wordFreq, nt);
    cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";
    
    return 0;
}