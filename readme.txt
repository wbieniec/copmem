copMEM, v0.2
============
copMEM is a program for efficient computation of MEMs (Maximal Exact Matches) in a pair of genomes.
Its main algorithmic idea requires that two internal parameters (k1 and k2) are coprime, hence the name.


Installation:
-------------
1) Linux
Use a 64-bit OS (e.g., Ubuntu / g++ 7.x).
Just "make all" in the directory with extracted files.
If you see any error messages, please contact the authors (wbieniec@kis.p.lodz.pl or sgrabow@kis.p.lodz.pl).

2) Windows (64-bit).
We recommend using Visual Studio 2017. The package contains the project file (.sln).
Of course, compile in Release 64-bit mode.


Usage:
------
./copmem -o <MEMs_file> [options] <R> <Q>

where the obligatory parameters are:
  -o <MEMs_file> output file with the found MEM information,
  <R>            reference file with one or more FASTA sequences,
  <Q>            query file with one or more FASTA sequences,

and the optional one:
  -l n           minimal MEM length (the default value is 100). This value must be >= 50,
  -v | -q        verbose mode (prints more details) | quiet mode,
  -f             fast mode, loads two files in memory, faster I/O operations,
  -r             calculates reverse complement matches. Does not go with -f switch,
  -b             calculates forward and reverse complement matches. Does not go with -f switch,
  -K 36|44|56    determines length of k-mer for hash calculation. Default is 44,
  -e             forces k2 == 1,
  -H 1|2|3|4|5   selects hash function. Default is 1.
     1: maRushPrime1HashSimplified (http://www.amsoftware.narod.ru/algo2.html);
     2: xxHash32 by Yann Collet (https://github.com/Cyan4973/xxHash);
     3: xxHash64 by Yann Collet (https://github.com/Cyan4973/xxHash);
     4: MetroHash64 by J. Andrew Rogers (https://github.com/jandrewrogers/MetroHash);
     5: CityHash64 from Google (https://github.com/google/cityhash).

copMEM prints maximal exact matches to the -o file.

Currently, copMEM is a stand-alone program. In a future version we are going to allow to use it also as a drop-in replacement for MUMmer.


MemoryFill:
-----------
MemoryFill is an auxiliary program that enables "cold start", especially for several runs with the same input files.
Usage:

./memoryfill <amount_of_memory>

Use with <amount_of_memory> close to the amount of the installed physical memory. One can use values succeeded by M or G, like 14G.
M = 2^20, G = 2^30.


Demo:
-----
Download genome archives: 
Homo sapiens ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
Mus musculus ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz

Unzip files into separate directories and concatenate single sequence files info one file, e.g.:
In Windows:
type mus\*.fa > mus.all.fa
type hum\*.fa > hum.all.fa
In Linux:
cat mus/*.fa > mus.all.fa
cat hum/*.fa > hum.all.fa

Run demo.cmd or demo.sh depending on the operating system.


Notes:
------
In the current version, copMEM is single-threaded. Despite that, it still works faster than its contenders (e.g., about twice faster than E-MEM using 10 threads on a 6c/12t i7 machine).
The memory usage of copMEM is, roughly, 
|R| + |maxQSeq| + 4 * (2^29 + |R| / k1) + |output|
bytes, where
maxQSeq is the longest sequence (i.e., chromosome) from file Q, |output| is the maximum output size *per sequence from file Q*, and k1 is an internal parameter approximately equal to sqrt(L - 44).  For example, k1 = 8 for L = 100 and k1 = 13 for L = 200.
If the datasets are larger than >4 GB (e.g., Triticum plants), replace the multiplier 4 with 8 in the memory formula above.
The value of 2^29 is the number of slots in the hash table.


Cite:
-----
If you find copMEM useful, please cite (at this moment, only an arXiv report exists):

https://arxiv.org/abs/1805.08816

Szymon Grabowski, Wojciech Bieniecki
"copMEM: Finding maximal exact matches via sampling both genomes"


References:
-----------
[E-MEM, an important competitor]
N. Khiste, L. Ilie 
E-MEM: efficient computation of maximal exact matches for very large genomes
Bioinformatics 31(4), 2015, pp. 509--514
(https://academic.oup.com/bioinformatics/article/31/4/509/2748225)

[essaMEM]
M. Vyverman, B. De Baets, V. Fack, P. Dawyndt
essaMEM: finding maximal exact matches using enhanced sparse suffix arrays
Bioinformatics 29(6), 2013, pp. 802--804
(https://academic.oup.com/bioinformatics/article/29/6/802/184020)
