# Vot
## A C++ library for computing variational optimal transportation

This package includes the prototype codes for reproducing the mapping in Figure 4 (b) of the paper:

Mi, Liang, Wen Zhang, Xianfeng Gu, and Yalin Wang. "Variational Wasserstein Clustering." In Proceedings of the European Conference on Computer Vision (ECCV), pp. 322-337. 2018.

![alt text](demo/sample.png?raw=true "Demo of Variational Wasserstein Clustering")

## Build and run

[![Build status](https://ci.appveyor.com/api/projects/status/7yw0ao44kfvjavfw?svg=true)](https://ci.appveyor.com/project/icemiliang/vot)

In the root directory, run:
```
$ make clean
$ make [-j4]
```

Check Makefile and config.mk to change paths if necessary.

A command-line program 'vot' is built from vot.cpp in demo/

A minimum command for running vot:
```
$ ./vot -e /data/empirical.vot -d /data/dirac.vot 
```

## Dependences
Boost 1.58, Eigen 3 (included in include/Eigen 3.3.5)

The code has been tested on Ubuntu 16.04 with g++ 5.5.0

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
  -h [ --help ]                         Print help messages
```

One iteration of the Wasserstein clustering is optimal transport. If you only want to compute the transport, set iterP to 1.

## Input and output files
Sample input files are provided in /demo/data/

### Input:

  *.vot specifies Dirac or empirical measures

  This program mainly deals with unstructured 3D point clouds. For 2D clouds, simply set 'z' values to zeros. It still regards the data as of 3D. This is due to the nature of Voro++ which is 3D.

  All samples must be within the range of (-1,1) in each dimension. The current version only supports a cubic boundary. Spherical boundary will be added soon.

  A Python version and an n-dimensional version is comming soon.

### Output:

  *.vot specifies the resulting Dirac measures.

  *.gnu specifies the resulting power Voronoi diagram.

Output files are compatible with Gnuplot. A sample Gnuplot script (gnuplotScripts.txt) is also provided in /demo/. If you are using the sample files, the script should plot the picture above.

## Code structure
The command-line program is built from /demo/vot.cpp

The main body of the Vot code is in /src/ot.cc and /src/diagram.cc

A static library is built in /lib/libvot.a

We give credit to Voro++ (http://math.lbl.gov/voro++/) for computing power Voronoi diagrams. Many files in /src/ are from Voro++, with our modification.

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
