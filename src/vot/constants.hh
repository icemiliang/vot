// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Oct 16th 2018

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

namespace vot {
    const double otTOLERANCE = 0.0000001;
    const double otTOTAL_MASS = 1.0;

    struct otBBox{
        const double xMin = -1.0 - otTOLERANCE;
        const double xMax = 1.0 + otTOLERANCE;
        const double yMin = -1.0 - otTOLERANCE;
        const double yMax = 1.0 + otTOLERANCE;
        const double zMin = -1.0 - otTOLERANCE;
        const double zMax = 1.0 + otTOLERANCE;
    }; // minX maxX minY maxY minZ maxZ

    const double otMaxDouble = std::numeric_limits<double>::max();
    const double otMinDouble = -std::numeric_limits<double>::max();

    const int METHOD_GD = 1; // Gradient descent
    const int METHOD_NEWTON = 2; // Newton's

    const double ballCenterX = 0, ballCenterY = 0, ballCenterZ = 0, ballRadius = 1.05;

    const int SUCCESS = 0;
    const int ERROR_IN_COMMAND_LINE = 1;
    const int ERROR_UNHANDLED_EXCEPTION = 2;

    const int ERR_INDEX_OUT_OF_RANGE = 11;
    const int ERR_DIMENSIONS_MISMATCH = 15;
    const int ERR_ZERO_DIVIDER = 21;

    const std::string ERR_MSG_INDEX_OUT_OF_RANGE = "[Error] Index out of range.";
    const std::string ERR_MSG_DIMENSION_MISMATCH = "[Error] Dimensions don't match.";
    const std::string ERR_MSG_ZERO_DIVIDER = "[Error] Divider too close to zero.";
}

#endif
