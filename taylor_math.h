
// By Ulf Ekstrom March-June 2009.
// Elementary functions. For inclusion in taylor.h only!

//Taylor series of 1/(a+x)
template<int N>
void inv_taylor(taylor<1,N>& t, const tnum_t &a)
{
  assert(a != tnum_t(0) && "1/(a+x) not analytic at a = 0");
  t[0] = 1/a;
  for (int i=1;i<=N;i++)
    t[i] = -t[i-1]/a;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<1,Ndeg> tmp;
  inv_taylor(tmp,t[0]);
  tmp*=x;
  taylor<Nvar, Ndeg> res;
  t.compose(res,tmp);
  return res;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp *= 1/x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const taylor<Nvar,Ndeg>&t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<1,Ndeg> tmp;
  inv_taylor(tmp,t2[0]);
  taylor<Nvar, Ndeg> res;
  t2.compose(res,tmp);
  res*=t1;
  return res;
}

// Evaluate the taylor series of exp(x0+x)=exp(x0)*exp(x)
template<int Ndeg>
void exp_taylor(taylor<1,Ndeg> &t, const tnum_t &x0)
{
  tnum_t ifac = 1;
  t[0] = exp(x0);
  for (int i=1;i<=Ndeg;i++)
    {
      ifac *= i;
      t[i] = t[0]/ifac;
    }
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> exp(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  exp_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

// Log series log(a+x) = log(1+x/a) + log(a)
template<int N>
void log_taylor(taylor<1,N> &t, const tnum_t &x0)
{
  assert(x0 != tnum_t(0) && "log(x) not analytic at x = 0");
  t[0] = log(x0);
  tnum_t xn = x0;
  for (int i=1;i<=N;i++)
    {
      t[i] = (2*(i & 1)-1)/(i*xn);
      xn *= x0;
    }
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> log(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  log_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template<int N>
void pow_taylor(taylor<1,N>& t, const tnum_t &x0, const tnum_t &a)
{
  assert(x0 != tnum_t(0) && "pow(x,a) not analytic at x = 0");
  t[0] = pow(x0,a);
  for (int i=1;i<=N;i++)
    t[i] = t[i-1]*(a-i+1)/(x0*i);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> pow(const taylor<Nvar,Ndeg> &t, const tnum_t &a)
{
  taylor<1,Ndeg> tmp;
  pow_taylor(tmp,t[0],a);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> sqrt(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  pow_taylor(tmp,t[0],0.5);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

//Integer exponent version is analytical at t[0] = 0
// TODO: speed this up a bit by using compose
template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> pow(const taylor<Nvar,Ndeg> &t, int n)
{
  if (n < 1)
    return pow(t,tnum_t(n));
  taylor<Nvar,Ndeg> res = t;
  while (n-- > 1)
      res *= t;
  return res;
}

// Use that d/dx atan(x) = 1/(1 + x^2),
// Taylor expand in x^2 and integrate.
template<int Ndeg>
void atan_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  // Calculate taylor expansion of 1/(1+a^2+x)
  taylor<1,Ndeg> invt,x;
  inv_taylor(invt,1+a*a);
  //insert x = 2*a*x + x^2
  x = taylor<1,Ndeg>(0);
  if (Ndeg > 0)
    x[1] = 2*a;
  if (Ndeg > 1)
    x[2] = 1;
  x.compose(t,invt);
  //Integrate each term and set the constant
  t.integrate();
  t[0] = atan(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> atan(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  atan_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


/* Taylor expansion of exp(-(a+x)^2) =
   exp(-a^2-2a*x)*exp(-x^2)

 */
template<int Ndeg>
void gauss_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  exp_taylor(t,-a*a);
  t.stretch(-2*a);
  taylor<1,Ndeg> g;
  g = 1;
  for (int i=1;i<=Ndeg/2;i++)
    g[2*i] = -g[2*(i-1)]/i;
  t*=g;
}

// Use that d/dx erf(x) = 2/sqrt(pi)*exp(-x^2),
// Taylor expand in x^2 and integrate.
template<int Ndeg>
void erf_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  gauss_taylor(t,a);
  t*=2/sqrt(TAYLOR_PI);
  //Integrate each term and set the constant
  t.integrate();
  t[0] = erf(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> erf(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  erf_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


template<int Ndeg>
void sin_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
    if (Ndeg > 0)
    {
	tnum_t s = sin(a), c = cos(a), fac = 1;
	for (int i=0;2*i<Ndeg;i++)
	{
	    t[2*i] = fac*s;
	    fac /= (2*i+1);
	    t[2*i+1] = fac*c;
	    fac /= -(2*i+2);
	}
	if (Ndeg % 2 == 0)
	    t[Ndeg] = s*fac;
    }
    else
    {
	t[0] = sin(a);
    }   
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> sin(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  sin_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


template<int Ndeg>
void cos_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
    if (Ndeg > 0)
    {
	tnum_t s = sin(a), c = cos(a), fac = 1;
	for (int i=0;2*i<Ndeg;i++)
	{
	    t[2*i] = fac*c;
	    fac /= -(2*i+1);
	    t[2*i+1] = fac*s;
	    fac /= (2*i+2);
	}
	if (Ndeg % 2 == 0)
	    t[Ndeg] = c*fac;
    }
    else
    {
	t[0] = cos(a);
    }   
}


template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> cos(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  cos_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


// hyperbolic arcsin function. d/dx asinh(x) = 1/sqrt(1+x^2)
// 1 + (a+x)^2 = 1+a^2 + 2ax + x^2
template<int Ndeg>
void asinh_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  taylor<1,Ndeg> tmp(1+a*a);
  if (Ndeg>0)
    tmp[1] = 2*a;
  if (Ndeg>1)
    tmp[2] = 1;
  t = pow(tmp,-0.5);
  t.integrate();
  t[0] = asinh(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> asinh(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  asinh_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


template<int Ndeg>
void sinc_taylor_at0(taylor<1,Ndeg>& t)
{
    t[0] = 1;
    tnum_t fac = 1;
    for (int i=1;i<=Ndeg/2;i++)
    {
	fac /= (2*i);
	t[2*i-1] = 0;
	fac /= -(2*i+1);
	t[2*i] = fac;
    }
    if (Ndeg % 2 == 1)
	t[Ndeg] = 0;
}


template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> sinc(const taylor<Nvar,Ndeg> &t)
{
    if (fabs(t[0]) < 1e-3)
    {
	taylor<1,Ndeg+8> tmp;
	sinc_taylor_at0(tmp);
	taylor<1,Ndeg> tmp2;	
	tnum_t dx = t[0];
	tmp.shift(tmp2,&dx);
	taylor<Nvar,Ndeg> res;
	t.compose(res,tmp2);
	return res;
    }
    else
    {
	return sin(t)/t;
    }
}


