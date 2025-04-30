# SSm0
Numerical implementation of the double-soft integral for massive-massless emitters at arbitrary angle

## Package structure
```
/
 |_ mmaSSm0                - Wolfram Library Link interface
 |_ tests
     |_ tests/IijEp0.dat   - quarks ep^1  term reference numbers
     |_ tests/IijEp-1.dat  - quarks ep^0  term reference numbers
     |_ tests/IijEp1.dat   - quarks ep^-1 term reference numbers
     |_ tests/IijEp-2.dat  - quarks ep^-2 term reference numbers
     |_ tests/SijEp0.dat   - gluons ep^1  term reference numbers
     |_ tests/SijEp-1.dat  - gluons ep^0  term reference numbers
     |_ tests/SijEp1.dat   - gluons ep^-1 term reference numbers
     |_ tests/SijEp-2.dat  - gluons ep^-2 term reference numbers
     |_ tests/test_SSm0.c  - Tests main program
 |_ ex_SSm0.c              - Example program
 |_ Makefile               
 |_ SSm0.h                 - Main library header file
 |_ SSm0_IijLiC.c          - Complex arg Li, quarks finite part
 |_ SSm0_IijLiR.c          - Real arg Li, quarks finite part
 |_ SSm0_internal.c        - Internal functions implementations
 |_ SSm0_internal.h        - Internal functions definitions
 |_ SSm0_li22.c            - Li22(x,y) function for real x,y
 |_ SSm0_Li22data.h        - Li22 log expansion tables
 |_ SSm0_Li22log1yNegTab.h - Li22(x,y) for x,y -> 1 and y < 0
 |_ SSm0_li2.c             - Real and complex argument Li2 function
 |_ SSm0_li3.c             - Real and complex argument Li3 function
 |_ SSm0_li4.c             - Real argument Li4 function
 |_ SSm0_SijLi22.c         - Li22 part of the gluons finite part
 |_ SSm0_SijLiC.c          - Complex arg Li, gluons finite part
 |_ SSm0_SijLiR.c          - Real arg Li, gluons finite part
 |_ SSm0_tld.c             - \tilde{S} and \tilde{I} functions
```

## Installation
A recent version of the package could be obtained from the git repo and compiled on any machine with a C compiler available
```
$ git clone https://github.com/apik/SSm0.git
$ cd SSm0
$ make
```
After compilation, the library file `libSSm0.a` and the example program `ex_SSm0` are produced

To check the validity of the code, we provide a set of reference points organized as tables in the format `beta cos(theta) value` for each expansion order,
to build and run tests, use
```
$ make check
$ cd tests
$ ./test_SSm0
```

## Usage
With the package, we provide an example program `ex_SSm0.c`, which calculates double soft integrals for the point specified as two command line arguments or in an arbitrary random point if there are no arguments
```
$ ./ex_SSm0 0.2 0.3
```
With produced output:
```
beta       =   0.20000000
cos(theta) =   0.30000000
                       \tilde{I}(quarks)                 \tilde{S}(gluons)
ep^-3                      0.0000000000                     -0.2500000000
ep^-2                     -0.0833333333                      0.0002022602
ep^-1                      0.1023516932                     -0.4409857886
ep^ 0                      0.0184460303                      0.0147284599
ep^ 1                      0.0847637908                      2.0817758455
```

To use compiled library `libSSm0.a` in user the user program like one present bellow
```C
/* useSSm0.c */
#include <stdio.h>
#include "SSm0.h"

int main()
{
  int epOrder = 1;
  double be   = 0.2;
  double ct   = 0.3;
  printf("%lf, %lf\n ", SSm0_tldS(epOrder, be, ct), SSm0_tldI(epOrder, be, ct));
  return 0;
}
```

one can compile and link it, assuming that the `SSm0` library can be found at `SSM0PATH` with
```
$ gcc -I${SSM0PATH} useSSm0.c -L${SSM0PATH} -lSSm0 -lm
```

## Wolfram Libray Link interface

In the folder `mmaSSm0`, we provide a simple interface to Mathematica, which can be built with

```
$ cd mmaSSm0
$ make
```
A simple test example for `x=0.2,y=0.3` can be run assuming that the Mathematica executable is `math` with
```
$ math -script test.m
```

More advanced example producing plots for `ep^1` parts of `\tilde{I}` and `\tilde{S}` can be run with
```
$ math -script plot.m
```
As a result, the two files `pltQQ.pdf` and `pltGG.pdf` should be produced


## References
This implementation is based on the calculation performed in 
 - Dennis Horstmann, Kirill Melnikov, Ming-Ming Long and Andrey Pikelner. Integral of the double-emission eikonal function for a massive and a massless emitter at an arbitrary angle. [arXiv:2504.20977](https://arxiv.org/abs/2504.20977)
