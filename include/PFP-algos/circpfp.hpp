/*
 * Multithread Prefix-free parse implementation to compute the non-circular PFP of sequence collections.
 * 
 * This code is adapted from https://github.com/alshai/Big-BWT/blob/master/newscan.hpp
 *
 */
#include <cstdint>
extern "C" {
#include "xerrors.h"
}
#include <vector>
#include <istream>
uint64_t kr_hash(string s);
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t terminator_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Res {
    size_t tot_char = 0;
    int us_th = 0;
    std::vector<uint64_t> terminator_positions;
};

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end; // input
  size_t parsed, words;  // output
  FILE *parse, *o;
  vector<uint64_t> terminator_positions; 
} mt_data;

// modified to handle terminators and non-circular parsing
void *non_circular_mt_parse_fasta(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d-> wordFreq;

  if(arg->verbose>1)
    printf("Scanning from %ld, size %ld as a FASTA record\n",d->true_start,d->true_end-d->true_start);
  if(d->true_end-d->true_start == 0) return NULL;

  // open input file
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }
  
  // prepare for parsing
  f.seekg(d->true_start); // move to the beginning of assigned region
  KR_window krw(arg->w);
  uint8_t c, pc = '\n';
  uint64_t base_pos = d->true_start; 
  uint64_t global_seq_number = 0; // Global sequence number counter
  
  // skip the header
  uint64_t current_pos = d->true_start;
  while((c = f.get()) != '\n' && current_pos < d->true_end) {
    current_pos++;
    // Count sequence headers we've passed to determine global sequence number
    if (c == '>')
      global_seq_number++;
  }
  pc = c;
  
  // parse the sequence(s)
  string current_word = "";
  bool found_first_trigger = false;
  uint64_t i = 0; // Position in the current sequence
  uint64_t last_trigger_pos = 0;
  uint64_t local_seq_number = 0; // Local sequence counter for this thread
  
  while((pc != EOF) && current_pos < d->true_end) {
    c = f.get();
    current_pos++;
    
    if(c > 64) { // A is 65 in ascii table
      c = std::toupper(c);
      if(c <= Dollar || c > 90) die("Invalid char found in input file. Exiting...");
      
      current_word.append(1, c);
      uint64_t hash = krw.addchar(c);
      
      // Check if this is a trigger position
      if (hash % arg->p == 0 && krw.current == arg->w) {
        if (!found_first_trigger) {
          // First trigger in this sequence
          found_first_trigger = true;
          
          // Save the position of the first trigger
          uint64_t start_char = i - arg->w + 1;
          if(fwrite(&start_char, sizeof(start_char), 1, d->o) != 1) 
            die("offset write error");
          
          // Add the first w-tuple to the dictionary
          uint64_t first_hash = kr_hash(current_word);
          if(fwrite(&first_hash, sizeof(first_hash), 1, d->parse) != 1) 
            die("parse write error");
          
          // Update frequency table
          xpthread_mutex_lock(&map_mutex, __LINE__, __FILE__);
          if (wordFreq->find(first_hash) == wordFreq->end()) {
            (*wordFreq)[first_hash].occ = 1;
            (*wordFreq)[first_hash].str = current_word;
          } else {
            (*wordFreq)[first_hash].occ += 1;
            if ((*wordFreq)[first_hash].str != current_word) {
              cerr << "Hash collision detected!\n";
              exit(1);
            }
          }
          xpthread_mutex_unlock(&map_mutex, __LINE__, __FILE__);
          
          d->words++;
          last_trigger_pos = i;
          current_word = current_word.substr(current_word.size() - arg->w);
        } else {
          // Not the first trigger, save the word from last trigger to here
          save_update_word(current_word, arg->w, *wordFreq, d->parse, false);
          d->words++;
          
          last_trigger_pos = i;
          current_word = current_word.substr(current_word.size() - arg->w);
        }
      }
      
      pc = c;
      i++;
    }
    else { 
      if (c == '>' || c == EOF || current_pos >= d->true_end) {
        // End of sequence, add terminator
        // Create a sequence-specific terminator with the global sequence number
        local_seq_number++;
        uint64_t seq_num = global_seq_number + local_seq_number;
        string terminator = "$";
        current_word.append(terminator);
        
        // Save terminator position
        uint64_t terminator_pos = base_pos + krw.tot_char;
        d->terminator_positions.push_back(terminator_pos);
        
        // Update the KR window with the terminator characters
        for (char tc : terminator) {
            krw.addchar(tc);
        }
        
        // Save the final word to the dictionary
        save_update_word(current_word, arg->w, *wordFreq, d->parse, true);
        d->words++;
        
        d->parsed += krw.tot_char;
        base_pos += krw.tot_char;
        
        // Reset for next sequence
        krw.reset();
        found_first_trigger = false;
        last_trigger_pos = 0;
        i = 0;
        current_word = "";
        
        // Skip header line
        while((c = f.get()) != EOF && current_pos < d->true_end) {
          current_pos++;
          if(c == '\n') break;
        }
        pc = c;
      }
    }
  }
  
  // If we ended with some partial sequence, finalize it
  if (current_word.size() > 0) {
    local_seq_number++;
    uint64_t seq_num = global_seq_number + local_seq_number;
    string terminator = "$";
    current_word.append(terminator);
    
    uint64_t terminator_pos = base_pos + krw.tot_char;
    d->terminator_positions.push_back(terminator_pos);
    
    for (char tc : terminator) {
        krw.addchar(tc);
    }
    
    save_update_word(current_word, arg->w, *wordFreq, d->parse, true);
    d->words++;
    d->parsed += krw.tot_char;
  }
  
  f.close(); 
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
Res parallel_parse_fasta(Args& arg, map<uint64_t,word_stats>& wf)
{
    assert(arg.th>0);
    pthread_t t[arg.th];
    mt_data td[arg.th];
    // scan file for start positions and execute threads
    FILE* fp = fopen(arg.inputFileName.c_str(), "r");
    if (fp == NULL) {
      throw new std::runtime_error("Cannot open input file " +arg.inputFileName);
    }
    fseek(fp, 0L, SEEK_END);
    size_t size = ftell(fp);
    rewind(fp);
    std::vector<size_t> th_sts(arg.th+1);
    th_sts[0] = 1; th_sts[arg.th] = size;
    for (int i = 1; i < arg.th; ++i) {
      th_sts[i] = (size_t) (size / arg.th) * i;
    }
  
    if(arg.verbose) {
    cout << "Thread: " << arg.th << endl;
    cout << "Total size: " << size << endl;
    cout << "------------------------" << endl; }
  
    uint8_t c = fgetc(fp); size_t hp = 1; bool sf = 0;
    size_t tstart = 1, tend = 1;
    size_t nt = 1;
    // this loop scans the Fasta file in order to properly divide it up
    // for the threads, so that they don't accidently start in a ">" header.
    // As soon as a proper start and end position has been found, execute the thread 
    for(int i=1;i<th_sts.size()-1;i++){
        
        size_t start = th_sts[i];
        size_t end = th_sts[i+1];
        fseek(fp, start, SEEK_SET);
        for(size_t j=start; j<end; ++j){
            c = fgetc(fp);
            if(c == '>'){hp=j;sf=1;break;}
        }
        
        if(sf==1){
            sf=0;
            //tend = hp-1;
            tend = hp;
            // prepare and execute thread j-1
            td[nt-1].wordFreq = &wf;
            td[nt-1].arg = &arg;
            td[nt-1].true_start = tstart;
            td[nt-1].true_end = tend;
            td[nt-1].words = 0;
            td[nt-1].parsed = 0;
            assert(td[nt-1].true_end <=size);
            td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
            td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
            xpthread_create(&t[nt-1],NULL,&non_circular_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);
            nt++;
            tstart = hp+1;
        }
    }
    
    tend = size;
    td[nt-1].wordFreq = &wf;
    td[nt-1].arg = &arg;
    td[nt-1].true_start = tstart;
    td[nt-1].true_end = tend;
    td[nt-1].words = 0;
    td[nt-1].parsed = 0;
    assert(td[nt-1].true_end <=size);
    td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
    td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
    xpthread_create(&t[nt-1],NULL,&non_circular_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);
    
    // wait for the threads to finish (in order) and close output files
    size_t tot_char=0;
    std::vector<uint64_t> all_terminator_positions;
    // Lock per proteggere l'accesso al vettore globale dei terminatori
    xpthread_mutex_lock(&terminator_mutex,__LINE__,__FILE__);
    
    for(int i=0;i<nt;i++) {
      xpthread_join(t[i],NULL,__LINE__,__FILE__); 
      all_terminator_positions.insert(all_terminator_positions.end(), 
      td[i].terminator_positions.begin(), 
      td[i].terminator_positions.end());
      if(arg.verbose) {
        cout << "s:" << td[i].true_start << "  e:" << td[i].true_end << "  pa:";
      }
      
      all_terminator_positions.insert(all_terminator_positions.end(), 
                              td[i].terminator_positions.begin(),
                              td[i].terminator_positions.end());
      
      // close thread-specific output files
      fclose(td[i].parse);
      fclose(td[i].o);
      if(td[i].words>0) {
        // extra check
        assert(td[i].parsed>arg.w);
        tot_char += td[i].parsed;
      }
      else assert(i>0); // the first thread must produce some words 
    }
    
    xpthread_mutex_unlock(&terminator_mutex,__LINE__,__FILE__);
    std::sort(all_terminator_positions.begin(), all_terminator_positions.end());
    
    fclose(fp);
    Res res; res.tot_char = tot_char; res.us_th = nt;
    res.terminator_positions = all_terminator_positions; 
    return res;
}