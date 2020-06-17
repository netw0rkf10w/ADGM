# Alternating Direction Graph Matching (ADGM)
This repository contains the code for ADGM introduced in the paper [Alternating Direction Graph Matching](https://khue.fr/papers/adgm_cvpr2017.pdf) (CVPR 2017) by [D. Khuê Lê-Huu](https://khue.fr) and [Nikos Paragios](http://cvn.centralesupelec.fr/~nikos).

A more recent implementation in Python can be found here: https://github.com/netw0rkf10w/pyADGM

v0.1, 24/02/2017, written by D. Khuê Lê-Huu.
 
If you use any part of this code, please cite:
```
@inproceedings{lehuu2017adgm,
 title={Alternating Direction Graph Matching},
 author={L{\^e}-Huu, D. Khu{\^e} and Paragios, Nikos},
 booktitle = {Proceedings of the {IEEE} Conference on Computer Vision and Pattern Recognition ({CVPR})},
 year = {2017}
}
```

## Note
- This is a preliminary re-implementation in C++ Eigen and should be 
considered as pre-release. I haven't tested it on the full benchmark yet 
and thus the performance is not guaranteed. If you observe some strange 
behavior then please let me know. 

- Furthermore, some caching part has not been re-implemented yet and thus, 
the current software might be slow (one can expect 2x-10x speedup in the 
future release version). If you compare running time in your paper, please 
indicate that.

- In the current version, only third-order potentials are supported.


## Installation
In Matlab, go to the folder ```ADGM/``` and run:
```
compile.m
```


## Usage
```
X = FUNCTION(X0, [], [], [], [], indH, valH, rho, MAX_ITER, verbose, eta, iter1, iter2);
```

where ```FUNCTION``` can be one of the following: ```ADGM1, ADGM2, ADGM1_SYMMETRIC, ADGM2_SYMMETRIC```.
If the third-order tensor valH is super-symmetric then you should use 
the ```_SYMMETRIC``` versions because they offer several times speedup by 
exploiting the symmetric structure of the tensor.

The parameters of the above function are:

Output X: the returned assignment matrix (N2 x N1)
Intput: 
- X0 the initial solution (N2 x N1)
- indH: matrix of dimension Nt x 3 representing the indices of the tensor valH (third-order)
- valH: vector of dimension Nt x 1 representing the values of the potential tensor
- rho: initial penalty parameter
- MAX_ITER: maximum number of iterations
- verbose: 'true' will print out the output
- eta (> 1.0), iter1, iter2 (denoted by eta, T1, T2 in the paper): parameters for 
the apdative scheme applied to the penalty parameter.

Typical values are:
rho = nP1*nP2/1000;
MAX_ITER = 5000;
eta = 2.0; 
iter1 = 200;
iter2 = 50;

+ Increasing rho or eta usually results in faster convergence but lower objective
(and vice-versa: decreasing them usually offer higher objective values) 
+ Decreasing iter1 or iter2 usually results in faster convergence but lower objective
(and vice-versa: increasing them usually offer higher objective values) 


**Important:**
- The indices in indH start from 0 (i.e. C++ indices and not Matlab ones)
- Although ADGM is formulated as a minimzation problem, the above function 
solves a MAXIMIZATION problem (for the ease of comparison with the other methods). 
Thus, the input potential tensor valH should represent the similarity between the graphs.


## Demo
The script demo.m implements a synthetic third-order graph matching problem 
and solves it using the two variants of ADGM as well as Duchenne's Tensor 
Matching algorithm (for comparison).

To succesfully run it, follow the steps below:

1. In Matlab, go to ```ann_mwrapper/``` and run:
```
ann_compile_mex
```

2. Go back to the main folder and run:
```
mex assignmentoptimal.cpp
```

3. Go to ```TM/``` and run:
```
mex GCC='/usr/bin/g++-4.8' mexSource/mexComputeFeature.cpp -output mex/mexComputeFeature
mex GCC='/usr/bin/g++-4.8' mexSource/mexTensorMatching.cpp -output mex/mexTensorMatching
```

4. Go back to the main folder and run:
```
demo
```

## Troubleshooting

When running the demo, if you obtain errors like this:
```
Invalid MEX-file '*.mexa64':
Missing symbol '_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_assignERKS4_' required by
'*.mexa64'
```
then probably your C++ compiler is not supported by MATLAB.

To resolve this issue,
you can run `mex -setup c++` and choose an appropritate compiler, as recommended by
[MATLAB's official documentation](https://fr.mathworks.com/help/matlab/matlab_external/choose-c-or-c-compilers.html).

Alternatively, you can specify a compiler as an argument for `mex`. For example, on Linux:
```
mex GCC='/usr/bin/g++-4.8' <file.cpp>
```
Make sure to replace `mex` with `mex GCC='/usr/bin/g++-4.8'` to all the commands in the previous sections (including those in `ADGM/compile.m`).

For any questions or bug reports, please send me an email.
