cHNLdecay
============

This is an implementation of partial widths described in the arXiv paper [arXiv:1805.08567v2](https://arxiv.org/abs/1805.08567v2) by Kyrylo Bondarenko, Alexey Boyarsky, Dmitry Gorbunov and Oleg Ruchayskiy. It uses the CERN scientific software toolkit [ROOT](https://root.cern.ch) and is written in C++ using the GNU Multiple Precision Arithmetic Library [GMP](https://gmplib.org).

## Installation
It is necessary to have GMP and ROOT installed on the system. Adapt the library paths in `Makefile` accordingly to point to the location of the libraries for GMP. ROOT should take care of this by itself as `root-config` is used. Additionally, if debugging is wished you need to install gperftools or uncomment the parts `-lprofiler -ltcmalloc` in `LDLIBS_DBG` in the `Makefile`.

Then simply run `make` in order to build the project. A binary called `cHNLdecay` is produced as a result.

## Usage
When running **cHNLdecay** you can parse various parameters to the program to configure the HNL.

The options are as follows:

```
--loglevel {debug,info,warning,error}
```
(one must be chosen, if omited `info` is the default)


```
--generations {11,13,15} 
```
(multiple can be passed separated by a `,` and no whitespace. The leptons are identified by their PDG ID) 

```
--angle FLOAT
```
(the value of the mixing angle `|U_alpha|^2`, NOTE the squared)

```
--mass INTEGER
```
(the mass of the HNL in MeV)


And an example would be:
```
./cHNLdecay --generations \\mu,e,\\tau --mass 20000 --angle 3.5e-6
```

## License
All rights reserved by **Fabian A.J. Thiele**, <fabian.thiele@cern.ch>
