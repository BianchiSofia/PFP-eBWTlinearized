# PFP-eBWT
This work presents the exemption of an existing implementation of the Prefix-Free Parsing technique for the computation of the Extended Burrows-Wheeler Transform (PFP-eBWT) in order to support non-circular genomic sequence collections, i.e. where the last character does not precede the first.

# Usage

### Construction of the eBWT:
```
usage: pfpebwt [-h] [-w WSIZE] [-p MOD] 
               input

Tool to compute the eBWT and the GCA of a string collection.

positional arguments:
  input                     input fasta file name

optional arguments:
  -h, --help                show this help message and exit
  -w WSIZE, --wsize WSIZE   sliding window size (def. 10)
  -p MOD, --mod MOD         hash modulus (def. 100)
```
The default PFP algorithm will run with one prime number and one remainder to search for the trigger strings. 
The `.info` file contains more information on the output files.

# Example
### Download and Compile

```console
git clone https://github.com/BianchiSofia/PFP-eBWTlinearized.git
cd PFP-eBWTlinearized
mkdir build
cd build
cmake ..
make
```

### Run on Example Data

```console
// Build the eBWT on a toy data set
python3 pfpebwtlinearized example.fasta 
// Build the eBWT with custom window and module
python3 pfpebwtlinearized -w 4 -p 10
```
# External resources

* [gSACA-K](https://github.com/felipelouza/gsa-is.git)
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
