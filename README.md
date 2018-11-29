# Vot
## A C++ library for computing variational optimal transportation

This package includes the prototype codes for reproducing the map in Figure 4 (b) of the paper:

Mi, Liang, Wen Zhang, Xianfeng Gu, and Yalin Wang. "Variational Wasserstein Clustering." In Proceedings of the European Conference on Computer Vision (ECCV), pp. 322-337. 2018.

![alt text](demo/sample.png?raw=true "Demo of Variational Wasserstein Clustering")

[Click me for a Python version](https://github.com/icemiliang/pyvot)

## Build and run

[![Build status](https://ci.appveyor.com/api/projects/status/7yw0ao44kfvjavfw?svg=true)](https://ci.appveyor.com/project/icemiliang/vot)

In the root directory, run:
```
$ make clean
$ make [-j4]
```

Check Makefile and config.mk to change paths if necessary.

We provide two versions of Vot -- Vot and Votx. Vot uses Newton's method with convex hulls and is for the 3-D Euclidean space. Votx uses gradient descent and is for n-D Euclidean space. If you want to get a 3D Voronoi diagram, use Vot. Otherwise, use Votx.

Two command-line program 'vot' and 'votx' are built from vot.cpp and votx.cpp in /demo/

A minimum command for running vot (votx):
```
$ ./vot -e /data/empirical.vot -d /data/dirac.vot
$ ./votx -e /data/empirical.votx -d /data/dirac.votx
```

## Dependences
Boost 1.58, Eigen 3 (included in include/Eigen 3.3.5)

The code has been tested on Ubuntu 16.04 with g++ 5.5.0 and macOS 10.13.6 with g++ 4.2.1 (Apple's LLVM 10.0.0).

Boost is mainly for I/O and is not necessary for using Vot. You can skip it if you want to use Vot as a library.

## Options
```
  -d [ --dirac ] arg                    Dirac measures
  -e [ --empirical ] arg                Empirical measures
  -p [ --iterP ] arg (=5)              Max iteration for updating positions of centroids/Dirac measures. 
                                        Default is 5.
  -h [ --iterH ] arg (=1000)            Max iteration for updating 'H', or size of cells. Default is 1000.
  -t [ --threshold ] arg (=1e-6)        Threshold for terminating the program. Default is 1e-06
  -s [ --scale ] arg (=1)               Scale of the output diagram. Default is 1.
  -r [ --rate ] arg (=0.2)              Learning rate for gradient descent. Default is 0.2. Not needed for
                                        Newtons's method.
  -o [ --outdir ] arg                   Output directory
     [ --help ]                         Print help messages
```

One iteration of the Wasserstein clustering is optimal transport. If you only want to compute the transport, set iterP to 1.

## Input and output files
Sample input files are provided in /demo/data/

### Input:
  *.vot and *.votx specifie Dirac or empirical measures

  If you want to use vot on 2D point clouds, simplying set the third coordinate to 0.

  All samples must be within the range of (-1,1) in each dimension. The total mass of Empirical measures must be equal to the total Dirac of Dirac measures.

### Output:

  *.vot or *votx specifies the resulting Dirac and Empirical measures.

  *.gnu specifies the resulting power Voronoi diagram from Vot.

  Output files of Vot are compatible with Gnuplot. A sample Gnuplot script (gnuplotScripts.txt) is also provided in /demo/. If you are using the sample files, the script should plot the picture above.

  Votx does not involve convex hulls. We provide a Matlab script to plot the resulging samples.

## Code structure
The command-line programs are built from /demo/vot.cpp and /demo/votx.cpp

The main body of the Vot code is in /src/vot/ot.cc and /src/vot/diagram.cc
For Votx. it is in /src/votx/otx.cc

A static library is built in /lib/libvot.a (/lib/libvotx.a)

## Reference
If you find the code helpful, please cite the following article:

Mi, Liang, Wen Zhang, Xianfeng Gu, and Yalin Wang. "Variational Wasserstein Clustering." In Proceedings of the European Conference on Computer Vision (ECCV), pp. 322-337. 2018.

```
@inproceedings{mi2018variational,
  title={Variational {W}asserstein Clustering},
  author={Mi, Liang and Zhang, Wen and Gu, Xianfeng and Wang, Yalin},
  booktitle={Proceedings of the European Conference on Computer Vision (ECCV)},
  pages={322--337},
  year={2018}
}
```

## Contact
Please contact Liang Mi icemiliang@gmail.com for any issues. 
