#ifndef TAYLOR_H
#define TAYLOR_H

#include "polymul.h"

#ifdef TAYLOR_QD
// Using the quad-double library QD.
// See http://crd.lbl.gov/~dhbailey/mpdist/
// Note that you _need_ to use fpu_fix_start/end on
// x86 machines for correct rounding.
#warning Using quad double precision
#include <qd/qd_real.h>
typedef qd_real tnum_t;
typedef qd_real tlongnum_t;
#define TAYLOR_PI qd_real::_pi
// QD does not provide an error function, so here is a slow and
// not so smart Taylor expansion.
static inline qd_real erf(const qd_real &x)
{
  qd_real sum = 0;
  qd_real x2 = x*x, term = x;
  int k = 0;  
  if (fabs(x) > 12) // erf(12) = 1 - 1e-64
    return 1;
  sum += term;
  if (fabs(x) > 1e-64)
    {
      double magn = erf(to_double(x)); // order of magnitude of answer.
      while (fabs(term/magn) > 1e-64)
	{
	  k++;
	  term *= qd_real(-(2*k-1)*x2)/(k*(2*k+1));
	  sum += term;
	}
    }
  return sum*2/sqrt(qd_real::_pi);
}
#else
// Standard numeric types. Note that long double is not always
// the same on all machines.
#include <cmath>
typedef double tnum_t; // Type of coefficients
typedef long double tlongnum_t; // Type used to accumulate sums
#define TAYLOR_PI M_PI
#endif

template<int Nvar, int Ndeg>
  class taylor : public polynomial<tnum_t, Nvar, Ndeg>
{
public:
  taylor(void) {}
  taylor(tnum_t c0) : polynomial<tnum_t, Nvar,Ndeg>(c0) {} 
  // Set the constant term to c0 and the first order term of
  // variable var to var_value.
  taylor(tnum_t c0, int var, tnum_t var_value = 1.0) : polynomial<tnum_t, Nvar,Ndeg>(c0) 
  {
    assert(var>=0);
    assert(var<Nvar);
    if (Ndeg>0)
      this->c[var+1] = var_value;
  }
  // Multiply each term with the correct factorials to get 
  // derivatives, i.e. x^2y^3 is multiplied by 2!3!
  void deriv_facs(void)
  {
    polydfac(*this);
  }
  void integrate(void)
  {
    assert(Nvar == 1 && "Implement this..");
    for (int i=Ndeg;i>0;i--)
      this->c[i] = this->c[i-1]/i;
    this->c[0] = 0.0;
  }
  // Set x := alpha*x, y = :alpha*y etc
  void stretch(tnum_t alpha)
  {
    tnum_t an = alpha;
    int k = 1;
    for (int i=1;i<=Ndeg;i++)
      {
	for (;k<polylen(Nvar,i);k++)
	  this->c[k] *= an;
	an *= alpha;
      }
  }
  void invert_parity(void)
  {
    this->stretch(-1); // TODO: optimized version of stretch(-1)
  }
  // Calculate a shifted version of this taylor expansion of one variable.
  // out(x) ~= this(x+dx). Of course the shifted taylor series is no longer
  // exact around the new "zero" point. For this reason NdegOut may be choosen
  // to be smaller than Ndeg.
  // TODO: implement for more variables, and without runtime polylen.
  template<int NdegOut>
  void shift(taylor<1,NdegOut> &out, const tnum_t dx[Nvar]) const
  {
      assert(Nvar == 1);
      assert(NdegOut <= Ndeg);
      tnum_t dxpow[Ndeg];
      out = 0;
      dxpow[0] = 1;
      for (int i=1;i<=Ndeg;i++)
	  dxpow[i] = dx[0]*dxpow[i-1];
      for (int i=0;i<=NdegOut;i++)
	  for (int j=i;j<=Ndeg;j++)
	      out[i] += polylen(j-i,i)*dxpow[j-i]*(*this)[j];
  }
  taylor<Nvar,Ndeg> operator-(void) const
    {
      taylor<Nvar,Ndeg> res = *this;
      for (int i=0;i<res.size;i++)
	res[i] *= -1;
      return res;
    }
#ifndef RVALUE_REF
  void operator-=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] -= t.c[i];
    }
  void operator+=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] += t.c[i];
    }
  void operator*=(const tnum_t& scale)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] *= scale;
    }
  template<int Ndeg2>
  void operator*=(const taylor<Nvar,Ndeg2>& t)
    {
      taylormul(*this,t);
    }
#else
  taylor<Nvar,Ndeg>& operator-=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] -= t.c[i];
      return *this;
    }
  taylor<Nvar,Ndeg>& operator+=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] += t.c[i];
      return *this;
    }
  taylor<Nvar,Ndeg>& operator*=(const tnum_t& scale)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] *= scale;
      return *this;
    }
  template<int Ndeg2>
  taylor<Nvar,Ndeg>& operator*=(const taylor<Nvar,Ndeg2>& t)
    {
      taylormul(*this,t);
      return *this;
    }
  void entrywise_mul(const taylor<Nvar, Ndeg> &t)
  {
    for (int i=0;i<this->size;i++)
      this->c[i] *= t[i];
  }
#endif
  /* Put sum_i coeff[i]*(this - this[0])^i in res,
     used when evaluating analytical functions of this
  */
  template<int Nres>
  void compose(taylor<Nvar,Nres>& res, const taylor<1,Nres> &coeff) const
    {
      assert(Nres >= Ndeg);
      taylor<Nvar,Ndeg> tmp = *this;
      tmp[0] = 0;
      res = 0;
      for (int i=Nres;i>0;i--)
	{
	  res[0] += coeff[i];
	  res *= tmp;
	}
      res[0] += coeff[0];
    }
  tnum_t dot(const taylor<Nvar, Ndeg> &t) const
  {
    tlongnum_t sum = 0;
    for (int i=0;i<this->size;i++)
      sum += this->c[i]*t.c[i];
    return sum;
  }
};

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp;
  taylormul(tmp,t1,t2);
  return tmp;
}



template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = t1[i]+t2[i];
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp = -t;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] -= x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = t1[i]-t2[i];
  return tmp;
}


#ifdef RVALUE_REF
template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator*(taylor<Nvar,Ndeg>&& t1, const taylor<Nvar,Ndeg>& t2)
{
  return t1*=t2;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator*(taylor<Nvar,Ndeg>&& t1, taylor<Nvar,Ndeg>&& t2)
{
  return t1*=t2;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator*(const taylor<Nvar,Ndeg>& t1, taylor<Nvar,Ndeg>&& t2)
{
  return t2*=t1;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator*(const tnum_t &x, taylor<Nvar,Ndeg>&& t)
{
  for (int i=0;i<t.size;i++)
    t[i] *= x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator*(taylor<Nvar,Ndeg>&& t, const tnum_t &x)
{
  for (int i=0;i<t.size;i++)
    t[i] *= x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> &&operator-(taylor<Nvar,Ndeg>&& t)
{
  for (int i=0;i<t.size;i++)
    t[i] *= -1;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator+(const tnum_t &x, taylor<Nvar,Ndeg>&& t)
{
  t[0] += x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator+(taylor<Nvar,Ndeg>&& t, const tnum_t &x)
{
  t[0] += x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator+(taylor<Nvar,Ndeg>&& t1, const taylor<Nvar,Ndeg>& t2)
{
  for (int i=0;i<t1.size;i++)
    t1[i] += t2[i];
  return t1;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator+(const taylor<Nvar,Ndeg>& t1, taylor<Nvar,Ndeg>&& t2)
{
  for (int i=0;i<t2.size;i++)
    t2[i] += t1[i];
  return t2;
}


template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator+(taylor<Nvar,Ndeg>&& t1, taylor<Nvar,Ndeg>&& t2)
{
  for (int i=0;i<t1.size;i++)
    t1[i] += t2[i];
  return t1;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator-(const tnum_t &x, taylor<Nvar,Ndeg>&& t)
{
  for (int i=0;i<t.size;i++)
    t[i] *= -1;
  t[0] += x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> &&operator-(taylor<Nvar,Ndeg>&& t, const tnum_t &x)
{
  t[0] -= x;
  return t;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator-(taylor<Nvar,Ndeg>&& t1, const taylor<Nvar,Ndeg>& t2)
{
  for (int i=0;i<t1.size;i++)
    t1[i] -= t2[i];
  return t1;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator-(const taylor<Nvar,Ndeg>& t1, taylor<Nvar,Ndeg>&& t2)
{
  for (int i=0;i<t2.size;i++)
    t2[i] = t1[i] - t2[i];
  return t2;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> && operator-(taylor<Nvar,Ndeg>&& t1, taylor<Nvar,Ndeg>&& t2)
{
  for (int i=0;i<t1.size;i++)
    t1[i] -= t2[i];
  return t1;
}

#endif



#include "taylor_math.h"

#endif
