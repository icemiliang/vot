// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef POINTBASE_HH
#define POINTBASE_HH

#include "vot.hh"

namespace vot {

class Point {
public:
    Point(double x, double y, double z) { mP[0] = x; mP[1] = y; mP[2] = z; };
    Point() { mP[0] = mP[1] = mP[2] = 0; };
    ~Point() {};

    double &operator[](int i)       { assert(0 <= i && i < 3); return mP[i]; };
    double  operator()(int i) const { assert(0 <= i && i < 3); return mP[i]; };
    double  operator[](int i) const { assert(0 <= i && i < 3); return mP[i]; };

    // l2-norm
    double norm() const { return sqrt(mP[0] * mP[0] + mP[1] * mP[1] + mP[2] * mP[2]); };
    double norm_squared() const { return mP[0] * mP[0] + mP[1] * mP[1] + mP[2] * mP[2]; };
    // lp-norm
    double norm(int i) const { return pow(pow(mP[0], i) + pow(mP[1], i) + pow(mP[2], i), 1 / double(i)); };

    Point &operator += (const Point &p) { mP[0] += p[0]; mP[1] += p[1]; mP[2] += p[2]; return *this; };
    Point &operator -= (const Point &p) { mP[0] -= p[0]; mP[1] -= p[1]; mP[2] -= p[2]; return *this; };
    Point &operator *= (const double s) { mP[0] *= s;    mP[1] *= s;    mP[2] *= s;    return *this; };
    Point &operator /= (const double s) { mP[0] /= s;    mP[1] /= s;    mP[2] /= s;    return *this; };

    Point operator+(const Point &p) const { Point np(mP[0] + p[0], mP[1] + p[1], mP[2] + p[2]); return np; };
    Point operator-(const Point &p) const { Point np(mP[0] - p[0], mP[1] - p[1], mP[2] - p[2]); return np; };
    Point operator*(const double s) const { Point np(mP[0] * s, mP[1] * s, mP[2] * s); return np; };
    Point operator/(const double s) const { Point np(mP[0] / s, mP[1] / s, mP[2] / s); return np; };
    Point operator-() const { Point p(-mP[0], -mP[1], -mP[2]); return p; };
    Point operator^(const Point & p) const {
        Point np(mP[1] * p[2] - mP[2] * p[1],
                 mP[2] * p[0] - mP[0] * p[2],
                 mP[0] * p[1] - mP[1] * p[0]);
        return np;
    };
    double operator *(const Point &p) const { return mP[0] * p[0] + mP[1] * p[1] + mP[2] * p[2]; };
    bool operator == (const Point &p) const { return (mP[0] == p[0] && mP[1] == p[1] && mP[2] == p[2]); };
    double angle(Point &p) { return acos((*this) * p / (norm() * p.norm())); };
    double x(){ return mP[0]; };
    double y(){ return mP[1]; };
    double z(){ return mP[2]; };

    friend std::ostream &operator<<(std::ostream &os, const Point pt) {
        os << pt[0] << " " << pt[1] << " " << pt[2];
        return os;
    }

    void cout() { std::cout << mP[0] << " " << mP[1] << " " << mP[2] << std::endl; };

protected:
    double mP[3];
};

}

#endif
