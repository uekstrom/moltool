#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>
#include "vec.h"

template<class numtype>
class quaternion : public vec<numtype,4>
{
protected:
  using vec<numtype,4>::c;
public:
  quaternion() {}
  quaternion(const numtype &c0) 
  {
    c[0] = c0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;
  }
  void operator=(const numtype &c0)
  {
    c[0] = c0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;
  }
  quaternion<numtype> conj(void) const
  {
    quaternion<numtype> q;
    q[0] = c[0];
    q[1] = -c[1];
    q[2] = -c[2];
    q[3] = -c[3];
    return q;
  }
  quaternion<numtype> operator*(const quaternion<numtype> &q) const
  {
    quaternion<numtype> r;
    r[0] = c[0]*q[0] - c[1]*q[1] - c[2]*q[2] - c[3]*q[3];
    r[1] = c[1]*q[0] + c[0]*q[1] - c[3]*q[2] + c[2]*q[3];
    r[2] = c[2]*q[0] + c[3]*q[1] + c[0]*q[2] - c[1]*q[3];
    r[3] = c[3]*q[0] - c[2]*q[1] + c[1]*q[2] + c[0]*q[3];
    return r;
  }
  void rotate(vec<numtype,3> &vout, const vec<numtype,3> &v) const
  {
    quaternion<numtype> qv;
    qv[0] = 0;
    for (int i=0;i<3;i++)
      qv[i+1] = v[i];
    qv = (*this)*qv*this->conj();
    for (int i=0;i<3;i++)
      vout[i] = qv[i+1];
  }
};

static inline double sinc(double x)
{
    if (fabs(x) > 1e-8)
	return sin(x)/x;
    else
	return 1 - x*x/6;
}

template<class numtype>
inline quaternion<numtype> exp(const quaternion<numtype> &q)
{
  quaternion<numtype> eq;
  numtype norm = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  eq[0] = cos(norm);
  numtype sn = sinc(norm);
  for (int i=1;i<4;i++)
    eq[i] = sn*q[i];
  eq *= exp(q[0]);
  return eq;
}

template<class numtype>
inline quaternion<numtype> log(const quaternion<numtype> &q)
{
  quaternion<numtype> eq;
  numtype normq = sqrt(q.abs2());
  numtype normv = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  eq[0] = log(normq);
  if (normv >  1e-12)
  {
      numtype sn = acos(q[0]/normq)/normv;
      for (int i=1;i<4;i++)
	  eq[i] = sn*q[i];
  }
  else
  {
      for (int i=1;i<4;i++)
	  eq[i] = 0;
  }
  return eq;
}

template<class numtype>
inline quaternion<numtype> quaternion_rotor(const vec<numtype,3> &v)
{
  quaternion<numtype> eq;
  numtype norm = sqrt(v.abs2());
  eq[0] = cos(norm);
  numtype sn = sinc(norm);
  for (int i=0;i<3;i++)
    eq[i+1] = sn*v[i];
  return eq;
}

#endif
