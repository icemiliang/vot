// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef OT_HH
#define OT_HH

#include "diagram.hh"

namespace vot {
    
class OT {
public:
    OT(){ otBBox bBox;  mDiagram = new Diagram(bBox, 3); };
    ~OT();

    void import_data(std::string pDiracFile, std::string pEmpiricalFile);
    // Import parameters
    void setup(const int pMaxIterP, const int pMaxIterH, const double pThres, 
               const double pLearnRate, const bool pVerbose, const double pPlotScale,
               const std::string pFilePrefix);
    void cluster();

protected:
    void read_empirical_from_file(std::string filename);
    void read_dirac_from_file(std::string filename);
    void check_total_mass();

    // Variables
    Diagram *mDiagram;

    std::string mOutFilePrefix;
    bool mFlagVerbose;
    int mMaxIterD;
    int mMaxIterH;
    double mThres;
    double mLearnRate;
    int mMethod; // Newton or GD
};

}
#endif
