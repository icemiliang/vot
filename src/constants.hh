// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef _constants_H_
#define _constants_H_

#include <limits>

namespace vot {

    const double otTOLERANCE = 0.00001;
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

    const int METHOD_GRADIENT = 1; // Gradient descent
    const int METHOD_NEWTON = 2; // Newton's

    const double ballCenterX = 0, ballCenterY = 0, ballCenterZ = 0, ballRadius = 1.05;

    const int SUCCESS = 0;
    const int ERROR_IN_COMMAND_LINE = 1;
    const int ERROR_UNHANDLED_EXCEPTION = 2;

    const int FUNCTION_VOT = 1;
    const int FUNCTION_VPM = 2; 
    const std::string FUNCTION_VOT_STRING = "Variational Optimal Transportation";
    const std::string FUNCTION_CLUSTER_STRING = "Variational Clustering";
    const std::string FUNCTION_SKELETON_STRING = "Variational Skeleton Extraction";
}

#endif