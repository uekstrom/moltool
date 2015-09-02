#ifndef TAYLOR_OPT_H
#define TAYLOR_OPT_H
#include "taylor.h"
#include "vec.h"


// Multiply the Hessian of the taylor expansion t by the vector x,
// put the result in y.
template<int Nvar>
static void hessvec(vec<tnum_t,Nvar> &y, const taylor<Nvar, 2> &t, const vec<tnum_t, Nvar> &x)
{
    int k = Nvar+1;
    for (int i=0;i<Nvar;i++)
    {
	y[i] = 0;
	for (int j=i;j<Nvar;j++)
	    y[i] += x[j]*t[k++];
    }
    k = Nvar + 1;
    for (int j=0;j<Nvar;j++)
	for (int i=j;i<Nvar;i++)
	    y[i] += x[j]*t[k++];
}

// Find the minimum x of a second order Taylor expansion. 
// The expansion must be positive definite to guarantee a solution.
// Return 0 if a solution was found, -1 if no solution was found, or 
// if a non-positive definite Hessian was detected.
// If no solution can be found consider using Levenberg regularization.
template<int Nvar>
int taylor_newton(vec<tnum_t,Nvar> &x, const taylor<Nvar, 2> &t, tnum_t thres)
{
    int n = 1;
    vec<tnum_t,Nvar> r, d, q;
    tnum_t delta_new, delta_old, stopping = -1, alpha;
    for (int i=0;i<Nvar;i++)
	x[i] = 0;
 restart:
    hessvec(r,t,x);
    if (r.dot(x) < 0) //Not positive definite
	return -1;
    for (int i=0;i<Nvar;i++)
	r[i] = -t[i+1] - r[i];
    d = r;
    delta_new = r.abs2();
    if (stopping < 0)
	stopping = delta_new*thres*thres;
    while (n < 2*Nvar and delta_new > stopping)
    {
	hessvec(q,t,d);
	tnum_t dq = d.dot(q);
	if (dq <= 0)
	    return -1;
	alpha = delta_new/dq;
	x += alpha*d;
	if (n % 50 != 0)
	    r -= alpha*q;
	else
	    goto restart;
	delta_old = delta_new;
	delta_new = r.abs2();
	d = r + (delta_new/delta_old)*d;
	n++;
    }
    if (delta_new > stopping)
	return -1;
    else
	return 0;
}


template<int Nvar>
int levenberg_step(vec<double, Nvar> &x, const taylor<Nvar,2> &t, double max_step)
{
    max_step *= max_step;
    double lambda = 0; // Levenberg shift
    int niter = 0;
    taylor<Nvar,2> devt = t;
    while (taylor_newton(x,devt,1e-14) != 0 or x.abs2() > max_step)
    {
	niter++;
	if (niter > 100)
	    return -1;
	// If Newton does not converge or the step is too large,
	// add a quadratic function around x = 0 to limit the step.
	// Remove old shift
	devt[4] -= lambda;
	devt[7] -= lambda;
	devt[9] -= lambda;
	if (niter == 1 and x.abs2() > max_step)
	{
	    // Calculate a lambda to make the step approximately
	    // of length max_step/sqrt(2).
	    double abs2g = 0;
	    for (int i=1;i<4;i++)
		abs2g += devt[i]*devt[i];
	    lambda = sqrt(abs2g/(2*max_step));
	}
	else if (niter == 1)
	    lambda = 1; // newton failed, probably because of non-definite hessian
	else
	    lambda *= 2;

	if (lambda > 1e10)
	    return -1;

	devt[4] += lambda;
	devt[7] += lambda;
	devt[9] += lambda;
    }
    return 0;
}


#endif
