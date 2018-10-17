// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Oct 16th 2018

#ifndef POINTBASEX_HH
#define POINTBASEX_HH

#include "basics.hh"

namespace votx {

// Point in the n-dimensioal vector space
class PointX {
public:
    //TODO: copy another PointX

    PointX(const int dim, double& coor) { vp = &coor; vn = dim;};

    PointX(const int dim) { vn = dim; vp = new double[vn]; std::fill(&vp[0],&vp[vn],0.0); }

    ~PointX() { delete [] vp; }; // destructor seems doing nothing.

    void print() const {
        std::cout << "(";
        for (int i = 0; i < vn-1; i++) {
            std::cout << vp[i] << ", ";
        }
        std::cout << vp[vn-1] << ")" << std::endl;
    };

    double operator[](int i) const { 
        ASSERT_VOT( i >= 0 && i < vn, ERR_MSG_INDEX_OUT_OF_RANGE); 
        return vp[i];
    };
    
    double &operator[](int i) { 
        ASSERT_VOT( i >= 0 && i < vn, ERR_MSG_INDEX_OUT_OF_RANGE); 
        return vp[i];
    };

    // Addition
    PointX operator + (const PointX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] + p[i]; }
        PointX res(vn,*tmp);
        return res;
    };
    // Subtraction
    PointX operator - (const PointX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] - p[i]; }
        PointX res(vn,*tmp);
        return res;
    };
    // Negative
    PointX operator - () const {
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = -vp[i]; }
        PointX res(vn,*tmp);
        return res;
    };
    // Element-wise product
    PointX operator * (const PointX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] * p[i]; }
        PointX res(vn,*tmp);
        return res;
    };
    // Scaling
    PointX operator * (const double s) const {
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] * s; }
        PointX res(vn,*tmp);
        return res;
    };
    // Scaling
    PointX operator / (const double s) const {
        ASSERT_VOT( fabs(s) > OTX_TOLERANCE, ERR_MSG_ZERO_DIVIDER);
        double* tmp = new double[vn];
        for (int i = 0; i < vn; i++) { tmp[i] = vp[i] / s; }
        PointX res(vn,*tmp);
        return res;
    };

    // In-place addition
    PointX &operator += (const PointX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] += p[i]; }
        return *this;
    };
    // in-pace substraction
    PointX &operator -= (const PointX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] -= p[i]; }
        return *this;
    };
    // In-place element-wise product
    PointX &operator *= (const PointX &p) {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        for (int i = 0; i < vn; i++) { vp[i] *= p[i]; }
        return *this;
    };
    // In-place scaling
    PointX &operator *= (const double s) {
        for (int i = 0; i < vn; i++) { vp[i] *= s; }
        return *this;
    };
    // In-place scaling
    PointX &operator /= (const double s) {
        ASSERT_VOT( fabs(s) > OTX_TOLERANCE, ERR_MSG_ZERO_DIVIDER);
        for (int i = 0; i < vn; i++) { vp[i] /= s; }
        return *this;
    };

    // Inner product
    double operator ^ (const PointX &p) const {
        ASSERT_VOT( p.dim() == vn, ERR_MSG_DIMENSION_MISMATCH);
        double sum = 0.0;
        for (int i = 0; i < vn; i++) { sum += vp[i] * p[i]; }
        return sum;
    };
    // Equality
    bool operator == (const PointX &p) const {
        bool equal = true;
        for (int i = 0; i < vn; i++) { equal = equal && (abs(vp[i] - p[i])<OTX_TOLERANCE); }
        return equal;
    }

    // l2-norm
    double norm() const { return sqrt(norm_p_powered(2)); };

    // l2-norm squared
    double norm_squared() const { return norm_p_powered(2); };

    // lp-norm
    double norm_p(double p) const { 
        if (fabs(p) < OTX_TOLERANCE) return 1.0;
        return pow(norm_p_powered(p), 1.0/p); 
    };

    // lp-norm squared
    double norm_p_powered(double p) const {
        if (fabs(p) < OTX_TOLERANCE) return 1.0;
        double res = 0.0;
        for (int i = 0; i < vn; i++) { res += pow(vp[i],p); }
        return res; 
    };

    // dist
    double distL22(const PointX &p) const {
        double dist = 0.0;
        for (int i = 0; i < vn; i++) {
            dist += (vp[i]-p[i]) * (vp[i]-p[i]); 
        }
        return dist;
    }
    double distL2(const PointX &p) const { return sqrt(distL22(p)); }
    
    int dim() const { return vn;};

    double* coor() { return vp; }

protected:

    double* vp; // pointer to point-array
    int vn;     // dimension
};

}

#endif
