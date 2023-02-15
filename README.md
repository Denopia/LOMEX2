## SAKE (Strobemer Assisted K-mer Extraction)

SAKE is a program for k-mer extraction. It utilizes two other programs: BFCounter (https://github.com/pmelsted/BFCounter) and SPOA (https://github.com/rvaser/spoa). Both programs are modified to fit the SAKE pipeline, so the original versions cannot be used. The modified versions are included in this repository. 

Small warning: SAKE is not optimized for performance or user-friendliness. SAKE is likely going to be slow with larger data sets. Due to it being a proof-of-concept prototype program, it is somewhat clunky to use.  


### How to install

1. Download the project directory.

2. Go to the main directory (SAKE).

3. Install the different programs of the pipeline in the following order
 
* To install modified BFCounter, run:

```
cd BFCounterForStrobemers
make
cd ..
```

* To install SPOA, run:

```
cd SAKEandSPOA/spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../../..
```

* To install SAKE, run:

```
cd SAKEandSPOA/sake
make sake
cd ../..
```
### How to use

* Input file only in fastq format.

1. Run Modified BFCounter count to find strobemers

```
./BFCounterForStrobemers/BFCounter count -n <int> -m <int> -l <int> -k <int> -w <int> -f <int> -t <int> -o <string> reads.fastq

Where:
n = approximate number of output strobemers (check the original BFCounter documentation)
m = number of minimizers
l = minimizer length
k = strobemer length (m*l)
w = minimizer window length
f = distance between two consecutive minimizer windows
t = number of threads
o = file path for the output strobemers (binary format)
reads.fastq is just the input read set file (fastq format)
```


2. Run Modified BFCounter dump to transform binary format strobemers into readable format
```
./BFCounterForStrobemers/BFCounter dump -k <int> -i <string> -o <string>

Where:
k = strobemer length
i = input binary file path (output from the previous step)
o = output file path for the readable strobemers
```

3. Run SAKE to produce consensus sequences
```
./SAKEandSPOA/sake/sake -k <string> -r <string> -o <string> -n <int> -l <int> -s <int> -f <int> -b <float> -z <float> -a <int> -g <int> -c <int> -t <int>

Where:
k = path to the file with readable strobemers (output from previous step)
r = path to the read set file (fastq format)
o = path to the output file where consensus sequences are written
n = number of minimizers in the strobemers
l = length of the minimizers in the strobemers
s = minimizer window size
f = distance between two consecutive minimizer windows
b = bundling threshold
z = relative bundle support threshold
a = minimum strobemer abundance
g = hard bundle support threshold
c = character support threshold
t = number of threads

```

4. Split consensus sequences into k-mers

```
python ./pytools/split_into_kmers.py <string1> <string2> <int>

Where:
The first string is the path to the file with the consensus sequences
The second string is the path to the output file where the k-mers are written
And the integer is the length of the k-mers
```

5. Orient k-mers into their canonical representation (lexicographic order)
```
python ./pytools/kmer_orienter.py <string1> <string2> <int1> <int2>

Where:
The first string is the path to the file with the k-mers
The second string is the path to the file where oriented k-mers are written
The first integer is the strobemer length
The second integer is an abundance threshold that should be set to 1 (This parameter does not make much sense anymore, it is due to leftover code from previous version of the script. Should be fixed at some point in the future.)

```

6. Sort and remove duplicate k-mers
```
sort <string1> > <string2> 
uniq <string2> > <string3>

Where:
String 1 is path to the file with oriented k-mers.
String 2 is path to the file where sorted and oriented k-mers are written.
String 3 is path to the file where unique sorted and oriented k-mers are written. [This is the final output file.]

```


### Dependencies

SAKE
* g++ 11.3.0 (older versions might also work)

SPOA
* gcc 4.8+ | clang 3.5+
* cmake 3.12+    
* (spoa_exe)(spoa_test) zlib 1.2.8+

SPOA Optional
* (optional) USCiLab/cereal 1.3.0
* (optional) simd-everywhere/simde 0.7.0
* (optional) google/cpu_features 0.6.0
* (spoa_exe)(spoa_test) rvaser/bioparser 3.0.13
* (spoa_exe)(spoa_test) rvaser/biosoup 0.10.0
* (spoa_test) google/googletest 1.10.0

### License
* SAKE is released under MIT license 
* Modified version of SPOA is released under its original license (MIT)
* Modified BFCounter is released under its original license (GPL3)
* (Also check the licenses of different parts of BFCounter in the original README file)
