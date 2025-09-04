/*
 * Code to build the eBWT and GCA of a text using its Prefix-free parse data structures.
 */

 #include <iostream>
 #include <string>
 #include <vector>
 #include <fstream>
 #include <ctype.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <ctime>
 #include <assert.h>
 #include <errno.h>
 #include <regex>
 #include <map>
 
 #include <sdsl/bit_vectors.hpp>
 #include <sdsl/int_vector.hpp>
 
 extern "C" {
 #include "utils.h"
 #include "xerrors.h" 
 }
 #include "common.hpp"
 #include "dictionary.hpp"
 #include "pfp_parse-gca.hpp"
 #include "pfp_parse.hpp"
 #include "pfp.hpp"
 #include "pfp_ebwt_gca.hpp" 
 #include "pfp_ebwt_gca_complete.hpp" 
 
 using namespace std;
 
 // struct containing command line parameters and other globals
 typedef struct {
    string inputFileName = "";
    int w = 10;
    bool rle = 0;
    bool complete = 0;
    bool sample = 0;
    bool sample_first = 0;
 } Args;
 
// Structure to represent an element of the modified Suffix Array
 struct ModifiedSAEntry {
     size_t original_index;  // Index in the original SA
     size_t position;        // Position in the dictionary
     std::string suffix;     // Corresponding suffix
     int terminal_number;    // Terminator number
 };
 
 static void parseArgs(int argc, char** argv, Args *arg ) {
   extern int optind, opterr, optopt;
   extern char *optarg;  
   int c;
 
   puts("==== Command line:"); 
   for(int i=0;i<argc;i++) 
     printf(" %s",argv[i]); 
   puts("\n");
 
   while ((c = getopt( argc, argv, "w:rsfc") ) != -1) { 
     switch(c) { 
       case 'w':
       arg->w = atoi(optarg); break; 
       case 'r':
       arg->rle = 1; break; 
       case 's':
       arg->sample = true; break;
       case 'f':
       arg->sample_first = true; break;
       case 'c':
       arg->complete = true; break;
       case '?':
       puts("Unknown option. Use -h for help.");
       exit(1);
     }
   }
 
   arg->inputFileName = argv[optind];
 }
 
 std::string extractCleanSuffix(const std::vector<uint8_t>& dict, size_t position, size_t max_length = 100) {
     std::string suffix = "";
     
     for (size_t j = position; j < std::min(position + max_length, dict.size()); j++) {
         if (dict[j] == EndOfWord || dict[j] == EndOfDict) {
             break;
         } else {
             suffix += static_cast<char>(dict[j]);
         }
     }
     
     return suffix;
 }
 
 int main(int argc, char** argv) {
     
     // translate command line parameters 
     Args arg;
     parseArgs(argc, argv, &arg);
     
     // start measuring wall time clock
     time_t start_wc = time(NULL);
 
     cout << "Computing BWT of the dictionary..." << endl;
     dictionary dict(arg.inputFileName,arg.w);
     
     cout << "Building the BWT of the dictionary took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
     start_wc = time(NULL);
 
     if( !arg.sample and !arg.complete )
     { 
       cout << "Loading parse's data structures..." << endl;
       pfp_parse pars(arg.inputFileName);
             
       cout << "Loading parse's data structures took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
       start_wc = time(NULL);
       
       // compute only the eBWT
       cout << "Computing eBWT of the text..." << endl;
       pfp pfp(pars,dict,arg.inputFileName,arg.w,arg.rle);
       cout << "Building the eBWT of Text took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
     
    }
     else if( arg.sample )
     {
       cout << "Loading parse's data structures..." << endl;
       pfp_parse_gca pars(arg.inputFileName);
       cout << "Loading parse's data structures took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
       start_wc = time(NULL);
     
       // check to verify terminators
       cout << "Verifying parse data structures for sequence terminators..." << endl;
       
       // compute the eBWT + gCA samples
       cout << "Computing the eBWT of the text and the GCA-samples..." << endl;
       try {
         pfp_ebwt_gca pfp_ebwt_gca(pars,dict,arg.inputFileName,arg.w,arg.rle,arg.sample_first);
         cout << "Building the eBWT and the GCA-samples of text took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
       } catch (const std::exception& e) {
         cerr << "Exception during eBWT/GCA construction: " << e.what() << endl;
         return 1;
       }
     }
     else if (arg.complete)
     {
       cout << "Loading parse's data structures..." << endl;
     
       try {
           pfp_parse_gca pars(arg.inputFileName);
           
           // Using pfp_ebwt_gca_complete with parameters that include terminators
           cout << "Computing the eBWT and the GCA of the text..." << endl;
           pfp_ebwt_gca_complete pfp_ebwt_gca_complete(pars, dict, arg.inputFileName, arg.w, arg.rle);
           
       } catch (const std::exception& e) {
           cerr << "Exception caught: " << e.what() << endl;
           return 1;
       } catch (...) {
           cerr << "Unknown exception caught!" << endl;
           return 1;
       }
     }
     
     return 0;
 }
