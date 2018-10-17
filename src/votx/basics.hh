// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Oct 16th 2018

#ifndef _BASICS_H_
#define _BASICS_H_

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
#include <limits>
#include <unordered_set>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace votx {
    #define ASSERT_VOT(condition, message) \
        do { \
            if (! (condition)) { \
                std::cerr << "--> Assertion `" #condition "` failed in " << __FILE__ \
                          << " line " << __LINE__ << ": " << message << std::endl; \
                std::terminate(); \
            } \
        } while (false);

    #define add_check_point(idx) \
        std::cout << " check point: " << idx << std::endl;

    const int SUCCESS = 0;
    const int ERROR_IN_COMMAND_LINE = 1;
    const int ERROR_UNHANDLED_EXCEPTION = 2;

    const int ERR_INDEX_OUT_OF_RANGE = 11;
    const int ERR_DIMENSIONS_MISMATCH = 15;
    const int ERR_ZERO_DIVIDER = 21;

    const std::string ERR_MSG_INDEX_OUT_OF_RANGE = "[Error] Index out of range.";
    const std::string ERR_MSG_DIMENSION_MISMATCH = "[Error] Dimensions don't match.";
    const std::string ERR_MSG_ZERO_DIVIDER = "[Error] Divider too close to zero.";

    const double OTX_TOLERANCE = 0.00000001;
    const double OTX_TOTAL_MASS = 1.0;

    const double OTX_MAX_DOUBLE = std::numeric_limits<double>::max();
    const double OTX_MIN_DOUBLE = -std::numeric_limits<double>::max();
}

#endif
