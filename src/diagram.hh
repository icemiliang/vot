// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef DIAGRAM_HH
#define DIAGRAM_HH

#include "voro++.hh"

#include "point.hh"

namespace vot {

class Diagram {
public:

    // Constructor requires a bounding box and # of blocks. If not provided, then [-1,1]^3.
    // # of blocks controls the complexity of finding the nearest particles to cut plane.
    Diagram(const otBBox bBox, const int numBlocks);
  
    // index x y z dirac mass h
    void add_dirac( double x, double y, double z, double dirac, double mass, double h, bool fix);

    // index x y z mass centroidIndex
    void add_empirical( double x, double y, double z, double m, int p);

    // Main process
    // Return success or error
    bool update(int method, double thres, double step, int iterD, int iterH); // update h
    
    // Return true if Dirac measures are stale
    bool update_dirac();
    bool update_dirac_beta();
    // Return max length between adjacent Dirac
    double interpolate();

    void write_results(int const iterD, int const iterH, std::string filePrefix);

    double total_dirac_mass();
    double total_empirical_mass();

    void setup(const bool verbose);

protected:
    void draw_diagram(std::string filename);
    void write_dirac(std::string filename);
    void write_dirac_as_input(std::string filename);
    void update_h(double step);
    void update_h(Eigen::SparseMatrix<double>& pHessian);

    double compute_wasserstein();
    double compute_gradient(); // Return the norm of gradients

    voro::container_poly *mDiagram;
    std::vector<Empirical> mEmpiricals;
    std::vector<Dirac> mDiracs;

    int mNumDiracs;
    int mNumEmpiricals;

    Eigen::VectorXd mGradient;

    std::vector<double> mHs;

    bool mFlagVerbose;
};

class DiagramSphere : public Diagram {
public:
    DiagramSphere(const otBBox bBox, const int numBlocks);
};

class DiagramCylinder : public Diagram {
public:
    DiagramCylinder(const otBBox bBox, const int numBlocks);
};
}

#endif