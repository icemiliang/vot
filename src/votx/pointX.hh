// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Oct 13th 2018

#ifndef POINTX_HH
#define POINTX_HH

#include "pointbaseX.hh"

namespace votx {

class EmpiricalX;
class DiracX;

// n-dimensional Dirac sample
class DiracX: public PointX {
public:
    DiracX(const int pDim, double& pCoor, const double pDirac, const double pMass)
    : PointX(pDim, pCoor) { vp = &pCoor; vn = pDim; vd = pDirac; vm = pMass; }

    DiracX(const int pDim, double& pCoor, const double pDirac)
    : PointX(pDim, pCoor) { vp = &pCoor; vn = pDim; vd = pDirac; vm = 0.0; } 

    DiracX(const int pDim, const int pDirac) 
    :PointX(pDim) {vn = pDim; vd = pDirac; vm = 0.0; 
                   vp = new double[vn]; std::fill(&vp[0],&vp[vn],0.0); }

    // algebra
    DiracX operator + (const DiracX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] + p[i]; }
        DiracX res(vn, *tmp, vd+p.dirac(), vm+p.mass());
        return res;
    }
    DiracX operator - (const DiracX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] - p[i]; }
        DiracX res(vn, *tmp, vd-p.dirac(), vm-p.mass());
        return res;
    }
    DiracX operator * (const double s) const {
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] * s; }
        DiracX res(vn, *tmp, vd, vm);
        return res;
    }
    DiracX operator / (const double s) const {
        ASSERT_VOT( fabs(s) > OTX_TOLERANCE, ERR_MSG_ZERO_DIVIDER);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] / s; }
        DiracX res(vn, *tmp, vd, vm);
        return res;
    }

    // in-place algebra
    DiracX &operator += (const DiracX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] += p[i]; }
        vm += p.mass(); vd += p.dirac();
        return *this;
    }
    DiracX &operator -= (const DiracX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] -= p[i]; }
        vm -= p.mass(); vd -= p.dirac();
        return *this;
    }
    DiracX &operator *= (const double s) {
        for (int i = 0; i < vn; i++) { vp[i] *= s; }
        return *this;
    }
    DiracX &operator /= (const double s) {
        ASSERT_VOT( fabs(s) > OTX_TOLERANCE, ERR_MSG_ZERO_DIVIDER);
        for (int i = 0; i < vn; i++) { vp[i] /= s; }
        return *this;
    }
    
    double dirac() const { return vd; }
    void set_dirac(const double pDirac) { vd = pDirac; }

    bool fix() const { return vf; }
    void set_fix(const bool pFix) { vf = pFix; }

    void add_empirical_idx(const int idx) { vIdxSet.insert(idx); }
    void clear_idx() { vIdxSet.clear(); };
    void reset_idx(std::unordered_set<int>& pSet) { vIdxSet.clear(); vIdxSet = pSet; };
    int num_empirical() const { return vIdxSet.size(); };

    void update_coor(double & pCoor) { delete [] vp; vp = &pCoor; }

    // void reset_coor() { std::fill(&vp[0],&vp[vn],0.0); }

    void print() {
        std::cout << "("  << "update? " << vic <<  "; dim" << vn << "; dirac: " << vd 
                  << "; mass: " << vm << "; coor: ";
        for (int i = 0; i < vn-1; i++) { std::cout << vp[i] << ", "; }
        std::cout << vp[vn-1] << ")" << std::endl;
    }

    double mass() const { return vm; }
    void set_mass(const double pMass) { vm = pMass; }
    void add_mass(const double pMass) { vm += pMass; }

    std::unordered_set<int> &set() { return vIdxSet; }

    bool update_flag() const { return vic && !vIdxSet.empty(); }
    void set_to_update() { vic = true; }
    void reset_upate_flag() { vic = false; }

protected:
    double vd; // dirac
    double vm;  // mass
    bool vf; // fixed coor
    bool vic; // empirical sample index changed

    std::unordered_set<int> vIdxSet; // Empirical idx
};

// n-dimensional Empirical sample
class EmpiricalX: public PointX {
public: 
    EmpiricalX(const int pDim, double& pCoor, const double pMass)
    :PointX(pDim, pCoor) { vp = &pCoor; vn = pDim; vm = pMass; vdi = -1; }

    EmpiricalX(const int pDim, const double pMass)
    :PointX(pDim) { vn = pDim; vm = pMass; vp = new double[vn]; std::fill(&vp[0],&vp[vn],0.0); vdi = -1; }

    EmpiricalX(const int pDim)
    :PointX(pDim) { vn = pDim; vm = 0.0; vp = new double[vn]; std::fill(&vp[0],&vp[vn],0.0); vdi = -1; }

    // double distL22(const DiracX &p) {
    //     ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
    //     double dist = 0;
    //     for (int i = 0; i < vn; i++) { dist += (vp[i]-p[i]) * (vp[i]-p[i]); }
    //     return dist;
    // }
    // double distL2(const DiracX &p) { return sqrt(distL22(p)); } 

    // false means idx doesn't change
    void set_dirac_idx(const int idx){ vdi = idx; }
    int di() const { return vdi; }

    void print() const {
        std::cout << "( dirac_idx: " << vdi << "; mass: " << vm << "; coor: ";
        for (int i = 0; i < vn-1; i++) {
            std::cout << vp[i] << ", ";
        }
        std::cout << vp[vn-1] << ")" << std::endl;
    }

    EmpiricalX &operator += (const EmpiricalX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] += p[i]; }
        vm += p.mass();
        return *this;
    }

    EmpiricalX operator * (const double s) const {
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] * s; }
        EmpiricalX res(vn,*tmp, vm);
        return res;
    }

    double mass() const { return vm; }

protected:
    int vdi; // dirac index
    double vm;  // mass
};

}

#endif
