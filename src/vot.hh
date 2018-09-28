// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef _VOT_H_
#define _VOT_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdint.h>
#include <cassert>
#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <iterator>

#include "constants.hh"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifndef NDEBUG
#   define ASSERT_VOT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "--> Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT_VOT(condition, message) do { } while (false)
#endif

#endif
