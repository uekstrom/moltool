#ifndef VEC_H
#define VEC_H

#include <cassert>

template<class numtype, int N>
class vec
{
 public:
  numtype c[N];

  vec() {}
  vec(const numtype &c0)
  {
    for (int i=0;i<N;i++)
      c[i] = c0;
  }
  template<class T2>
  vec(const vec<T2,N> &src)
  {
      for (int i=0;i<N;i++)
	  c[i] = src.c[i];
  }
  vec<numtype,N> &operator=(const numtype &c0)
  {
    for (int i=0;i<N;i++)
      c[i] = c0;
    return *this;
  }
  const numtype &operator[](int i) const
  {
    assert(i>=0);
    assert(i<N);
    return c[i];
  }
  numtype &operator[](int i)
  {
    assert(i>=0);
    assert(i<N);
    return c[i];
  }
  // NB: multiplied with c using *=, watch out if * not commutative
  // for numtype.
  void operator*=(const numtype c0)
  {
    for (int i=0;i<N;i++)
      c[i] *= c0;
  }
  void operator+=(const numtype c0)
  {
    for (int i=0;i<N;i++)
      c[i] += c0;
  }
  void operator+=(const vec<numtype,N> &v)
  {
    for (int i=0;i<N;i++)
      c[i] += v[i];
  }
  void operator-=(const numtype c0)
  {
    for (int i=0;i<N;i++)
      c[i] -= c0;
  }
  void operator-=(const vec<numtype,N> &v)
  {
    for (int i=0;i<N;i++)
      c[i] -= v[i];
  }
  vec<numtype,N> operator+(const vec<numtype,N> &v) const
  {
    vec<numtype,N> r;
    for (int i=0;i<N;i++)
      r[i] = c[i] + v[i];
    return r;
  }
  vec<numtype,N> operator-(const vec<numtype,N> &v) const
  {
    vec<numtype,N> r;
    for (int i=0;i<N;i++)
      r[i] = c[i] - v[i];
    return r;
  }
  vec<numtype,N> operator-(void) const
  {
    vec<numtype,N> r;
    for (int i=0;i<N;i++)
      r[i] = -c[i];
    return r;
  }
  numtype abs2(void) const
  {
    numtype sum = 0;
    for (int i=0;i<N;i++)
      sum += c[i]*c[i];
    return sum;
  }
  vec<numtype, N> normalized(void) const
  {
    return 1/sqrt(abs2())**this;
  }
  numtype dot(const vec<numtype,N> &v) const
  {
    numtype sum = 0;
    for (int i=0;i<N;i++)
      sum += c[i]*v[i];
    return sum;
  }
};

template<class numtype, int N>
vec<numtype,N> operator*(const numtype &c, const vec<numtype,N> &v)
{
  vec<numtype,N> r;
  for (int i=0;i<N;i++)
    r[i] = c*v[i];
  return r;
}

// We need a version for c from the left or the right because numtype
// may not commute under *.
template<class numtype, int N>
vec<numtype,N> operator*(const vec<numtype,N> &v, const numtype &c)
{
  vec<numtype,N> r;
  for (int i=0;i<N;i++)
    r[i] = v[i]*c;
  return r;
}



#endif
