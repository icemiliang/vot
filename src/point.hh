// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef POINT_HH
#define POINT_HH

#include "pointbase.hh"

namespace vot {

// Dirac samples
class Dirac: public Point {
public:
    Dirac(const double pX, const double pY, const double pZ, 
             const double pDirac, const double pMass) {
        mP[0] = pX;  mP[1] = pY;  mP[2] = pZ; mDirac = pDirac; mMass = pMass;
    };
    Dirac(const double pX, const double pY, const double pZ,
        const double pDirac, const double pMass, const bool fix) {
        mP[0] = pX;  mP[1] = pY;  mP[2] = pZ; mDirac = pDirac; mMass = pMass; mFix = fix;
    };

    void reset_dirac(const double pDirac) { mDirac = pDirac; };
    void reset_mass() { mMass = 0; };
    void set_mass(const double pMass) { mMass = pMass; };
    void add_mass(double massFromEmpirical) { mMass += massFromEmpirical; };
    void set_fix(const bool fix) { mFix = fix; };
    
    bool fix() const { return mFix; };
    double dirac() const { return mDirac; };
    double mass() const { return mMass; };
    Dirac &operator -= (const Dirac &p) {
        mP[0] -= p[0]; mP[1] -= p[1]; mP[2] -= p[2]; mMass -= p.mass(); 
        return *this; 
    };
    Dirac operator*(const double s) const {
        Dirac np(mP[0] * s, mP[1] * s, mP[2] * s, mDirac, mMass); 
        return np; 
    };
    Dirac &operator += (const Dirac &p) { 
        mP[0] += p[0]; mP[1] += p[1]; mP[2] += p[2]; mDirac += p.dirac(); mMass += p.mass(); 
        return *this; 
    };
    Dirac &operator /= (const double s) { 
        assert(fabs(s) > otTOLERANCE); mP[0] /= s; mP[1] /= s; mP[2] /= s; //mMass /= s; 
        return *this; 
    };
    // Have to redefine the overloaded operators?
    Dirac operator-(const Dirac &p) const { 
        Dirac np(mP[0] - p[0], mP[1] - p[1], mP[2] - p[2], 0.0, 0.0); 
        return np; 
    };

    void cout() { std::cout << mP[0] << " " << mP[1] << " " << mP[2] << " " << mMass << std::endl; };
    bool operator ==(const Dirac &p) const {
        for (int i = 0; i < 2; i++) {
            if (abs(mP[i] - p[i]) > otTOLERANCE) return false;
        }
        return true;
    };
protected:

    double mDirac;
    double mMass;
    bool mFix;
};


// Empirical samples
class Empirical : public Point {
public:
    // Constructor
    // Initialize from scratch
    Empirical(const double pX, const double pY, const double pZ, const double pMass, const int pCellIdx) {
         mP[0] = pX;  mP[1] = pY;  mP[2] = pZ; mMass = pMass; mCellIndex = pCellIdx;
    };

    const double mass() { return mMass; };
    const int cellIndex() { return mCellIndex; };
    void set_cellIndex(const int cellIdx) { mCellIndex = cellIdx; };
    void cout() { std::cout << mP[0] << " " << mP[1] << " " << mP[2] << " " << mMass << std::endl; };
    // Transportation cost from empirical to dirac
    double costL2(Dirac &pDirac) const {
       return mMass *((mP[0] - pDirac[0]) * (mP[0] - pDirac[0]) +
                      (mP[1] - pDirac[1]) * (mP[1] - pDirac[1]) +
                      (mP[2] - pDirac[2]) * (mP[2] - pDirac[2]));
    };

    double costL2_tmp(Dirac &pDirac) const {
       return ((mP[0] - pDirac[0]) * (mP[0] - pDirac[0]) +
               (mP[1] - pDirac[1]) * (mP[1] - pDirac[1]) +
               (mP[2] - pDirac[2]) * (mP[2] - pDirac[2]));
    };
    
protected:
    //double mP[3]; // Derived from Point
    double mMass;
    int mCellIndex;

};


}

#endif
