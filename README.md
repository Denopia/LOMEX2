## SAKE (Strobemer Assisted K-mer Extraction)

SAKE is a program for k-mer extraction. It utilizes two other programs: BFCounter (https://github.com/pmelsted/BFCounter) and SPOA (https://github.com/rvaser/spoa). Both programs are modified to fit the SAKE pipeline, so the original versions cannot be used. The modified versions are included in this repository. 


### How to install

1. Download the project directory.

2. Go to the main directory "SAKE".

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


### Dependencies

* Check the BFCounter and SPOA dependencies.
* (Maybe something more?)

### License
* SAKE is released under MIT license 
* Modified version of SPOA is released under its original license (MIT)
* Modified BFCounter is released under its original license (GPL3)
* (Also check the licenses of different parts of BFCounter in the original README file)
