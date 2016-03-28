#ifdef __IBXLINCOA_H__
#ifndef __IBXLINCOA_H__
#define __IBXLINCOA_H__

static Int32 c__1 = 1;

namespace nsOptimisation
{

template<typename tFunc,Int32 tRows> int LincoaSearch(tFunc &lpFunction,Int32 npt, Int32 m, Float64 
	*a, Int32 *ia, Float64 *b, cMatrix<Float64,tRows,1> &x, Float64 rhobeg, 
	Float64 rhoend, Int32 *iprint, Int32 maxfun)
{

    /* System generated locals */
    Int32 a_dim1, a_offset, i__1, i__2;
    Float64 d__1;

    /* Local variables */
    Int32 i__, j, ib, np, iw, iac, irc, igo, iqf, irf, ihq, ixb, ifv,ipq, isp, ixn, ixo, ixp, ixs;
    Float64 sum;
    Int32 ndim;
    Float64 temp;
    Int32 nptm;
    Float64 zero;
    Int32 istp, ipqw, iflag, iamat, ibmat, izmat;
    Float64 smallx;

	UInt32 wSize = M*(2+tRows + npt*(4+tRows+npt) + tRows*(9+3*tRows);
	UInt32 TwSize;
	TwSize = max(M+3*tRows, 2*M+tRows);
	wSize+= max(TwSize,2*npt);
	Float64 *w = [wSize];

/*     This subroutine seeks the least value of a function of many variables, */
/*       subject to general linear inequality constraints, by a trust region */
/*       method that forms quadratic models by interpolation. Usually there */
/*       is much freedom in each new model after satisfying the interpolation */
/*       conditions, which is taken up by minimizing the Frobenius norm of */
/*       the change to the second derivative matrix of the model. One new */
/*       function value is calculated on each iteration, usually at a point */
/*       where the current model predicts a reduction in the least value so */
/*       far of the objective function subject to the linear constraints. */
/*       Alternatively, a new vector of variables may be chosen to replace */
/*       an interpolation point that may be too far away for reliability, and */
/*       then the new point does not have to satisfy the linear constraints. */
/*       The arguments of the subroutine are as follows. */

/*     N must be set to the number of variables and must be at least two. */
/*     NPT must be set to the number of interpolation conditions, which is */
/*       required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices */
/*       of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be */
/*       highly inefficent when the number of variables is substantial, due */
/*       to the amount of work and extra difficulty of adjusting more points. */
/*     M must be set to the number of linear inequality constraints. */
/*     A is a matrix whose columns are the constraint gradients, which are */
/*       required to be nonzero. */
/*     IA is the first dimension of the array A, which must be at least N. */
/*     B is the vector of right hand sides of the constraints, the J-th */
/*       constraint being that the scalar product of A(.,J) with X(.) is at */
/*       most B(J). The initial vector X(.) is made feasible by increasing */
/*       the value of B(J) if necessary. */
/*     X is the vector of variables. Initial values of X(1),X(2),...,X(N) */
/*       must be supplied. If they do not satisfy the constraints, then B */
/*       is increased as mentioned above. X contains on return the variables */
/*       that have given the least calculated F subject to the constraints. */
/*     RHOBEG and RHOEND must be set to the initial and final values of a */
/*       trust region radius, so both must be positive with RHOEND<=RHOBEG. */
/*       Typically, RHOBEG should be about one tenth of the greatest expected */
/*       change to a variable, and RHOEND should indicate the accuracy that */
/*       is required in the final values of the variables. */
/*     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the */
/*       amount of printing. Specifically, there is no output if IPRINT=0 and */
/*       there is output only at the return if IPRINT=1. Otherwise, the best */
/*       feasible vector of variables so far and the corresponding value of */
/*       the objective function are printed whenever RHO is reduced, where */
/*       RHO is the current lower bound on the trust region radius. Further, */
/*       each new value of F with its variables are output if IPRINT=3. */
/*     MAXFUN must be set to an upper bound on the number of calls of CALFUN, */
/*       its value being at least NPT+1. */
/*     W is an array used for working space. Its length must be at least */
/*       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ]. */
/*       On return, W(1) is set to the final value of F, and W(2) is set to */
/*       the total number of function evaluations plus 0.5. */

/*     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set */
/*       F to the value of the objective function for the variables X(1), */
/*       X(2),...,X(N). The value of the argument F is positive when CALFUN */
/*       is called if and only if the current X satisfies the constraints */
/*       to working accuracy. */

/*     Check that N, NPT and MAXFUN are acceptable. */

    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --x;
    --w;

    /* Function Body */
    zero = 0.;
    smallx = *rhoend * 1e-6;
    np = *tRows + 1;
    nptm = *npt - np;
    if (*tRows <= 1) {
		GAI::LogError("(/4x,\002Return from LINCOA because N is less than 2.\002)");
	
	goto L80;
    }
    if (*npt < *tRows + 2 || *npt > (*tRows + 2) * np / 2) {
	GAI::LogError("Return from LINCOA because NPT is not in\002,\002 the required interval.");
	
	goto L80;
    }
    if (*maxfun <= *npt) {
	GAI::LogError("Return from LINCOA because MAXFUN is less\002,\002 than NPT+1");
	
	goto L80;
    }

/*     Normalize the constraints, and copy the resultant constraint matrix */
/*       and right hand sides into working space, after increasing the right */
/*       hand sides if necessary so that the starting point is feasible. */

/* Computing MAX */
    i__1 = *m + *tRows * 3, i__2 = (*m << 1) + *tRows, i__1 = max(i__1,i__2), i__2 = *
	    npt << 1;
    iamat = max(i__1,i__2) + 1;
    ib = iamat + *m * *tRows;
    iflag = 0;
    if (*m > 0) {
	iw = iamat - 1;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = zero;
	    temp = zero;
	    i__2 = *tRows;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum += a[i__ + j * a_dim1] * x[i__];
 /*L40: */
/* Computing 2nd power */
		d__1 = a[i__ + j * a_dim1];
		temp += d__1 * d__1;
	    }
	    if (temp == zero) {
		GAI::LogError("Return from LINCOA because the gradient of\002,\002 a constraint is zero.");
		
		goto L80;
	    }
	    temp = sqrt(temp);
	    if (sum - b[j] > smallx * temp) {
		iflag = 1;
	    }
/* Computing MAX */
	    d__1 = b[j];
	    w[ib + j - 1] = max(d__1,sum) / temp;
	    i__2 = *tRows;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++iw;
/* L60: */
		w[iw] = a[i__ + j * a_dim1] / temp;
	    }
	}
    }
    if (iflag == 1) {
	if (*iprint > 0) {
	    GAI::LogWarning("LINCOA has made the initial X feasible by\002,\002 increasing part(s) of B.");
	    
	}
    }

/*     Partition the working space array, so that different parts of it can be */
/*     treated separately by the subroutine that performs the main calculation. */

    ndim = *npt + *tRows;
    ixb = ib + *m;
    ixp = ixb + *tRows;
    ifv = ixp + *tRows * *npt;
    ixs = ifv + *npt;
    ixo = ixs + *tRows;
    igo = ixo + *tRows;
    ihq = igo + *tRows;
    ipq = ihq + *tRows * np / 2;
    ibmat = ipq + *npt;
    izmat = ibmat + ndim * *tRows;
    istp = izmat + *npt * nptm;
    isp = istp + *tRows;
    ixn = isp + *npt + *npt;
    iac = ixn + *tRows;
    irc = iac + *tRows;
    iqf = irc + *m;
    irf = iqf + *tRows * *tRows;
    ipqw = irf + *tRows * np / 2;

/*     The above settings provide a partition of W for subroutine LINCOB. */

    Linco_Do<tRows>(lpFunction,tRows, npt, m, &w[iamat], &w[ib], &x[1], rhobeg, rhoend, iprint, 
	    maxfun, &w[ixb], &w[ixp], &w[ifv], &w[ixs], &w[ixo], &w[igo], &w[
	    ihq], &w[ipq], &w[ibmat], &w[izmat], &ndim, &w[istp], &w[isp], &w[
	    ixn], &w[iac], &w[irc], &w[iqf], &w[irf], &w[ipqw], &w[1]);
L80:
		delete []w;
return 0;
} /* lincoa_ */


template<Int32 tRows,typename tFunc> int Linco_Do(tFunc &lpFunction,Int32 *n, Int32 *npt, Int32 *m, Float64 
									 *amat, Float64 *b, cMatrix<Float64,tRows,1> &x, Float64 rhobeg, Float64 
									 rhoend, Int32 *iprint, Int32 *maxfun, Float64 *xbase, 
									 Float64 *xpt, Float64 *fval, Float64 *xsav, Float64 *xopt,
									 Float64 *gopt, Float64 *hq, Float64 *pq, Float64 *bmat, 
									 Float64 *zmat, Int32 *ndim, Float64 *step, Float64 *sp, 
									 Float64 *xnew, Int32 *iact, Float64 *rescon, Float64 *qfac,
									 Float64 *rfac, Float64 *pqw, Float64 *w)
{

	/* System generated locals */
	Int32 amat_dim1, amat_offset, xpt_dim1, xpt_offset, bmat_dim1, 
		bmat_offset, zmat_dim1, zmat_offset, qfac_dim1, qfac_offset, i__1,
		i__2, i__3;
	Float64 d__1, d__2;


	/* Local variables */
	Float64 f;
	Int32 i__, j, k, ih, nf, nh, ip, np;
	Float64 del, one;
	Int32 idz;
	Float64 rho, sum, ssq, diff, half;
	Int32 nact, knew;
	Float64 temp, fopt;
	Int32 kopt, nptm;
	Float64 zero, sumz;
	Int32 ifeas;
	Float64 delta, xdiff;
	Int32 nvala, nvalb;
	Float64 fsave;
	Int32 ksave;
	Float64 ratio, vquad, tenth, vqalt;
	Int32 itest;
	Float64 snorm, dffalt;
	Float64 delsav;
	Float64 distsq;
	Float64 qoptsq, xoptsq;


	/*     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are */
	/*       identical to the corresponding arguments in SUBROUTINE LINCOA. */
	/*     AMAT is a matrix whose columns are the constraint gradients, scaled */
	/*       so that they have unit length. */
	/*     B contains on entry the right hand sides of the constraints, scaled */
	/*       as above, but later B is modified for variables relative to XBASE. */
	/*     XBASE holds a shift of origin that should reduce the contributions */
	/*       from rounding errors to values of the model and Lagrange functions. */
	/*     XPT contains the interpolation point coordinates relative to XBASE. */
	/*     FVAL holds the values of F at the interpolation points. */
	/*     XSAV holds the best feasible vector of variables so far, without any */
	/*       shift of origin. */
	/*     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of */
	/*       the feasible vector of variables that provides the least calculated */
	/*       F so far, this vector being the current trust region centre. */
	/*     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT. */
	/*     HQ holds the explicit second derivatives of the quadratic model. */
	/*     PQ contains the parameters of the implicit second derivatives of the */
	/*       quadratic model. */
	/*     BMAT holds the last N columns of the big inverse matrix H. */
	/*     ZMAT holds the factorization of the leading NPT by NPT submatrix */
	/*       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, */
	/*       where the elements of DZ are plus or minus one, as specified by IDZ. */
	/*     NDIM is the first dimension of BMAT and has the value NPT+N. */
	/*     STEP is employed for trial steps from XOPT. It is also used for working */
	/*       space when XBASE is shifted and in PRELIM. */
	/*     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT, */
	/*       followed by STEP^T XPT(K,.), K=1,2,...,NPT. */
	/*     XNEW is the displacement from XBASE of the vector of variables for */
	/*       the current calculation of F, except that SUBROUTINE TRSTEP uses it */
	/*       for working space. */
	/*     IACT is an Int32 array for the indices of the active constraints. */
	/*     RESCON holds useful information about the constraint residuals. Every */
	/*       nonnegative RESCON(J) is the residual of the J-th constraint at the */
	/*       current trust region centre. Otherwise, if RESCON(J) is negative, the */
	/*       J-th constraint holds as a strict inequality at the trust region */
	/*       centre, its residual being at least |RESCON(J)|; further, the value */
	/*       of |RESCON(J)| is at least the current trust region radius DELTA. */
	/*     QFAC is the orthogonal part of the QR factorization of the matrix of */
	/*       active constraint gradients, these gradients being ordered in */
	/*       accordance with IACT. When NACT is less than N, columns are added */
	/*       to QFAC to complete an N by N orthogonal matrix, which is important */
	/*       for keeping calculated steps sufficiently close to the boundaries */
	/*       of the active constraints. */
	/*     RFAC is the upper triangular part of this QR factorization, beginning */
	/*       with the first diagonal element, followed by the two elements in the */
	/*       upper triangular part of the second column and so on. */
	/*     PQW is used for working space, mainly for storing second derivative */
	/*       coefficients of quadratic functions. Its length is NPT+N. */
	/*     The array W is also used for working space. The required number of */
	/*       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA. */

	/*     Set some constants. */

	/* Parameter adjustments */
	qfac_dim1 = *n;
	qfac_offset = 1 + qfac_dim1;
	qfac -= qfac_offset;
	amat_dim1 = *n;
	amat_offset = 1 + amat_dim1;
	amat -= amat_offset;
	zmat_dim1 = *npt;
	zmat_offset = 1 + zmat_dim1;
	zmat -= zmat_offset;
	xpt_dim1 = *npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	--b;
	--x;
	--xbase;
	--fval;
	--xsav;
	--xopt;
	--gopt;
	--hq;
	--pq;
	bmat_dim1 = *ndim;
	bmat_offset = 1 + bmat_dim1;
	bmat -= bmat_offset;
	--step;
	--sp;
	--xnew;
	--iact;
	--rescon;
	--rfac;
	--pqw;
	--w;

	/* Function Body */
	half = .5;
	one = 1.;
	tenth = .1;
	zero = 0.;
	np = *n + 1;
	nh = *n * np / 2;
	nptm = *npt - np;

	/*     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT, */
	/*       ZMAT and SP for the first iteration. An important feature is that, */
	/*       if the interpolation point XPT(K,.) is not feasible, where K is any */
	/*       Int32 from [1,NPT], then a change is made to XPT(K,.) if necessary */
	/*       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT */
	/*       is set so that XPT(KOPT,.) is the initial trust region centre. */

	Linco_Initialise(lpFunction,n, npt, m, &amat[amat_offset], &b[1], &x[1], rhobeg, iprint, &
		xbase[1], &xpt[xpt_offset], &fval[1], &xsav[1], &xopt[1], &gopt[1]
	, &kopt, &hq[1], &pq[1], &bmat[bmat_offset], &zmat[zmat_offset], &
		idz, ndim, &sp[1], &rescon[1], &step[1], &pqw[1], &w[1]);

	/*     Begin the iterative procedure. */

	nf = *npt;
	fopt = fval[kopt];
	rho = *rhobeg;
	delta = rho;
	ifeas = 0;
	nact = 0;
	itest = 3;
L10:
	knew = 0;
	nvala = 0;
	nvalb = 0;

	/*     Shift XBASE if XOPT may be too far from XBASE. First make the changes */
	/*       to BMAT that do not depend on ZMAT. */

L20:
	fsave = fopt;
	xoptsq = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* L30: */
		/* Computing 2nd power */
		d__1 = xopt[i__];
		xoptsq += d__1 * d__1;
	}
	if (xoptsq >= delta * 1e4 * delta) {
		qoptsq = xoptsq * .25;
		i__1 = *npt;
		for (k = 1; k <= i__1; ++k) {
			sum = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L40: */
				sum += xpt[k + i__ * xpt_dim1] * xopt[i__];
			}
			sum -= half * xoptsq;
			w[*npt + k] = sum;
			sp[k] = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				xpt[k + i__ * xpt_dim1] -= half * xopt[i__];
				step[i__] = bmat[k + i__ * bmat_dim1];
				w[i__] = sum * xpt[k + i__ * xpt_dim1] + qoptsq * xopt[i__];
				ip = *npt + i__;
				i__3 = i__;
				for (j = 1; j <= i__3; ++j) {
					/* L50: */
					bmat[ip + j * bmat_dim1] = bmat[ip + j * bmat_dim1] + 
						step[i__] * w[j] + w[i__] * step[j];
				}
			}
		}

		/*     Then the revisions of BMAT that depend on ZMAT are calculated. */

		i__3 = nptm;
		for (k = 1; k <= i__3; ++k) {
			sumz = zero;
			i__2 = *npt;
			for (i__ = 1; i__ <= i__2; ++i__) {
				sumz += zmat[i__ + k * zmat_dim1];
				/* L60: */
				w[i__] = w[*npt + i__] * zmat[i__ + k * zmat_dim1];
			}
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
				sum = qoptsq * sumz * xopt[j];
				i__1 = *npt;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L70: */
					sum += w[i__] * xpt[i__ + j * xpt_dim1];
				}
				step[j] = sum;
				if (k < idz) {
					sum = -sum;
				}
				i__1 = *npt;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L80: */
					bmat[i__ + j * bmat_dim1] += sum * zmat[i__ + k * 
						zmat_dim1];
				}
			}
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				ip = i__ + *npt;
				temp = step[i__];
				if (k < idz) {
					temp = -temp;
				}
				i__2 = i__;
				for (j = 1; j <= i__2; ++j) {
					/* L90: */
					bmat[ip + j * bmat_dim1] += temp * step[j];
				}
			}
		}

		/*     Update the right hand sides of the constraints. */

		if (*m > 0) {
			i__2 = *m;
			for (j = 1; j <= i__2; ++j) {
				temp = zero;
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L100: */
					temp += amat[i__ + j * amat_dim1] * xopt[i__];
				}
				/* L110: */
				b[j] -= temp;
			}
		}

		/*     The following instructions complete the shift of XBASE, including the */
		/*       changes to the parameters of the quadratic model. */

		ih = 0;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			w[j] = zero;
			i__1 = *npt;
			for (k = 1; k <= i__1; ++k) {
				w[j] += pq[k] * xpt[k + j * xpt_dim1];
				/* L120: */
				xpt[k + j * xpt_dim1] -= half * xopt[j];
			}
			i__1 = j;
			for (i__ = 1; i__ <= i__1; ++i__) {
				++ih;
				hq[ih] = hq[ih] + w[i__] * xopt[j] + xopt[i__] * w[j];
				/* L130: */
				bmat[*npt + i__ + j * bmat_dim1] = bmat[*npt + j + i__ * 
					bmat_dim1];
			}
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			xbase[j] += xopt[j];
			xopt[j] = zero;
			/* L140: */
			xpt[kopt + j * xpt_dim1] = zero;
		}
	}

	/*     In the case KNEW=0, generate the next trust region step by calling */
	/*       TRSTEP, where SNORM is the current trust region radius initially. */
	/*       The final value of SNORM is the length of the calculated step, */
	/*       except that SNORM is zero on return if the projected gradient is */
	/*       unsuitable for starting the conjugate gradient iterations. */

	delsav = delta;
	ksave = knew;
	if (knew == 0) {
		snorm = delta;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L150: */
			xnew[i__] = gopt[i__];
		}
		trstep_(n, npt, m, &amat[amat_offset], &b[1], &xpt[xpt_offset], &hq[1]
		, &pq[1], &nact, &iact[1], &rescon[1], &qfac[qfac_offset], &
			rfac[1], &snorm, &step[1], &xnew[1], &w[1], &w[*m + 1], &pqw[
				1], &pqw[np], &w[*m + np]);

				/*     A trust region step is applied whenever its length, namely SNORM, is at */
				/*       least HALF*DELTA. It is also applied if its length is at least 0.1999 */
				/*       times DELTA and if a line search of TRSTEP has caused a change to the */
				/*       active set. Otherwise there is a branch below to label 530 or 560. */

				temp = half * delta;
				if (xnew[1] >= half) {
					temp = delta * .1999;
				}
				if (snorm <= temp) {
					delta = half * delta;
					if (delta <= rho * 1.4) {
						delta = rho;
					}
					++nvala;
					++nvalb;
					temp = snorm / rho;
					if (delsav > rho) {
						temp = one;
					}
					if (temp >= half) {
						nvala = (Int32) zero;
					}
					if (temp >= tenth) {
						nvalb = (Int32) zero;
					}
					if (delsav > rho) {
						goto L530;
					}
					if (nvala < 5 && nvalb < 3) {
						goto L530;
					}
					if (snorm > zero) {
						ksave = -1;
					}
					goto L560;
				}
				nvala = (Int32) zero;
				nvalb = (Int32) zero;

				/*     Alternatively, KNEW is positive. Then the model step is calculated */
				/*       within a trust region of radius DEL, after setting the gradient at */
				/*       XBASE and the second derivative parameters of the KNEW-th Lagrange */
				/*       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively. */

	} else {
		/* Computing MAX */
		d__1 = tenth * delta;
		del = max(d__1,rho);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L160: */
			w[i__] = bmat[knew + i__ * bmat_dim1];
		}
		i__1 = *npt;
		for (k = 1; k <= i__1; ++k) {
			/* L170: */
			pqw[k] = zero;
		}
		i__1 = nptm;
		for (j = 1; j <= i__1; ++j) {
			temp = zmat[knew + j * zmat_dim1];
			if (j < idz) {
				temp = -temp;
			}
			i__2 = *npt;
			for (k = 1; k <= i__2; ++k) {
				/* L180: */
				pqw[k] += temp * zmat[k + j * zmat_dim1];
			}
		}
		qmstep_(n, npt, m, &amat[amat_offset], &b[1], &xpt[xpt_offset], &xopt[
			1], &nact, &iact[1], &rescon[1], &qfac[qfac_offset], &kopt, &
				knew, &del, &step[1], &w[1], &pqw[1], &w[np], &w[np + *m], &
				ifeas);
	}

	/*     Set VQUAD to the change to the quadratic model when the move STEP is */
	/*       made from XOPT. If STEP is a trust region step, then VQUAD should be */
	/*       negative. If it is nonnegative due to rounding errors in this case, */
	/*       there is a branch to label 530 to try to improve the model. */

	vquad = zero;
	ih = 0;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
		vquad += step[j] * gopt[j];
		i__1 = j;
		for (i__ = 1; i__ <= i__1; ++i__) {
			++ih;
			temp = step[i__] * step[j];
			if (i__ == j) {
				temp = half * temp;
			}
			/* L190: */
			vquad += temp * hq[ih];
		}
	}
	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		temp = zero;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			temp += xpt[k + j * xpt_dim1] * step[j];
			/* L200: */
			sp[*npt + k] = temp;
		}
		/* L210: */
		vquad += half * pq[k] * temp * temp;
	}
	if (ksave == 0 && vquad >= zero) {
		goto L530;
	}

	/*     Calculate the next value of the objective function. The difference */
	/*       between the actual new value of F and the value predicted by the */
	/*       model is recorded in DIFF. */

L220:
	++nf;
	if (nf > *maxfun) {
		--nf;
		if (*iprint > 0) {
			GAI::LogError("Return from LINCOA because CALFUN has been\002,\002 called MAXFUN times.");

		}
		goto L600;
	}
	xdiff = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xnew[i__] = xopt[i__] + step[i__];
		x[i__] = xbase[i__] + xnew[i__];
		/* L240: */
		/* Computing 2nd power */
		d__1 = x[i__] - xsav[i__];
		xdiff += d__1 * d__1;
	}
	xdiff = sqrt(xdiff);
	if (ksave == -1) {
		xdiff = rho;
	}
	if (xdiff <= tenth * rho || xdiff >= delta + delta) {
		ifeas = 0;
		if (*iprint > 0) {
			GAI::LogError("Return from LINCOA because rounding errors\002,\002 prevent reasonable changes to X.");

		}
		goto L600;
	}
	if (ksave <= 0) {
		ifeas = 1;
	}
	f = (Float64) ifeas;
	lpFunction(n, &x[1], &f);
	if (*iprint == 3) {
		GAI::LogError("Function number\002,i6,\002    F =\002,1pd18.10,\002    The corresponding X is:\002/(2x,5d15.6))");
		do_fio(&c__1, (Char *)&nf, (ftnlen)sizeof(Int32));
		do_fio(&c__1, (Char *)&f, (ftnlen)sizeof(Float64));
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			do_fio(&c__1, (Char *)&x[i__], (ftnlen)sizeof(Float64));
		}

	}
	if (ksave == -1) {
		goto L600;
	}
	diff = f - fopt - vquad;

	/*     If X is feasible, then set DFFALT to the difference between the new */
	/*       value of F and the value predicted by the alternative model. */

	if (ifeas == 1 && itest < 3) {
		i__1 = *npt;
		for (k = 1; k <= i__1; ++k) {
			pqw[k] = zero;
			/* L270: */
			w[k] = fval[k] - fval[kopt];
		}
		i__1 = nptm;
		for (j = 1; j <= i__1; ++j) {
			sum = zero;
			i__2 = *npt;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L280: */
				sum += w[i__] * zmat[i__ + j * zmat_dim1];
			}
			if (j < idz) {
				sum = -sum;
			}
			i__2 = *npt;
			for (k = 1; k <= i__2; ++k) {
				/* L290: */
				pqw[k] += sum * zmat[k + j * zmat_dim1];
			}
		}
		vqalt = zero;
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			sum = zero;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				/* L300: */
				sum += bmat[k + j * bmat_dim1] * step[j];
			}
			vqalt += sum * w[k];
			/* L310: */
			vqalt += pqw[k] * sp[*npt + k] * (half * sp[*npt + k] + sp[k]);
		}
		dffalt = f - fopt - vqalt;
	}
	if (itest == 3) {
		dffalt = diff;
		itest = 0;
	}

	/*     Pick the next value of DELTA after a trust region step. */

	if (ksave == 0) {
		ratio = (f - fopt) / vquad;
		if (ratio <= tenth) {
			delta = half * delta;
		} else if (ratio <= .7) {
			/* Computing MAX */
			d__1 = half * delta;
			delta = max(d__1,snorm);
		} else {
			temp = sqrt(2.) * delta;
			/* Computing MAX */
			d__1 = half * delta, d__2 = snorm + snorm;
			delta = max(d__1,d__2);
			delta = min(delta,temp);
		}
		if (delta <= rho * 1.4) {
			delta = rho;
		}
	}

	/*     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point */
	/*       can be moved. If STEP is a trust region step, then KNEW is zero at */
	/*       present, but a positive value is picked by subroutine UPDATE. */

	update_(n, npt, &xpt[xpt_offset], &bmat[bmat_offset], &zmat[zmat_offset], 
		&idz, ndim, &sp[1], &step[1], &kopt, &knew, &pqw[1], &w[1]);
	if (knew == 0) {
		if (*iprint > 0) {
			GAI::LogError("Return from LINCOA because the denominator of the updating formula is zero.");

		}
		goto L600;
	}

	/*     If ITEST is increased to 3, then the next quadratic model is the */
	/*       one whose second derivative matrix is least subject to the new */
	/*       interpolation conditions. Otherwise the new model is constructed */
	/*       by the symmetric Broyden method in the usual way. */

	if (ifeas == 1) {
		++itest;
		if (abs(dffalt) >= tenth * abs(diff)) {
			itest = 0;
		}
	}

	/*     Update the second derivatives of the model by the symmetric Broyden */
	/*       method, using PQW for the second derivative parameters of the new */
	/*       KNEW-th Lagrange function. The contribution from the old parameter */
	/*       PQ(KNEW) is included in the second derivative matrix HQ. W is used */
	/*       later for the gradient of the new KNEW-th Lagrange function. */

	if (itest < 3) {
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			/* L330: */
			pqw[k] = zero;
		}
		i__2 = nptm;
		for (j = 1; j <= i__2; ++j) {
			temp = zmat[knew + j * zmat_dim1];
			if (temp != zero) {
				if (j < idz) {
					temp = -temp;
				}
				i__1 = *npt;
				for (k = 1; k <= i__1; ++k) {
					/* L340: */
					pqw[k] += temp * zmat[k + j * zmat_dim1];
				}
			}
			/* L350: */
		}
		ih = 0;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			w[i__] = bmat[knew + i__ * bmat_dim1];
			temp = pq[knew] * xpt[knew + i__ * xpt_dim1];
			i__1 = i__;
			for (j = 1; j <= i__1; ++j) {
				++ih;
				/* L360: */
				hq[ih] += temp * xpt[knew + j * xpt_dim1];
			}
		}
		pq[knew] = zero;
		i__1 = *npt;
		for (k = 1; k <= i__1; ++k) {
			/* L370: */
			pq[k] += diff * pqw[k];
		}
	}

	/*     Include the new interpolation point with the corresponding updates of */
	/*       SP. Also make the changes of the symmetric Broyden method to GOPT at */
	/*       the old XOPT if ITEST is less than 3. */

	fval[knew] = f;
	sp[knew] = sp[kopt] + sp[*npt + kopt];
	ssq = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xpt[knew + i__ * xpt_dim1] = xnew[i__];
		/* L380: */
		/* Computing 2nd power */
		d__1 = step[i__];
		ssq += d__1 * d__1;
	}
	sp[*npt + knew] = sp[*npt + kopt] + ssq;
	if (itest < 3) {
		i__1 = *npt;
		for (k = 1; k <= i__1; ++k) {
			temp = pqw[k] * sp[k];
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L390: */
				w[i__] += temp * xpt[k + i__ * xpt_dim1];
			}
		}
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L400: */
			gopt[i__] += diff * w[i__];
		}
	}

	/*     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the */
	/*       least calculated value so far with a feasible vector of variables. */

	if (f < fopt && ifeas == 1) {
		fopt = f;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			xsav[j] = x[j];
			/* L410: */
			xopt[j] = xnew[j];
		}
		kopt = knew;
		snorm = sqrt(ssq);
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
			if (rescon[j] >= delta + snorm) {
				rescon[j] = snorm - rescon[j];
			} else {
				rescon[j] += snorm;
				if (rescon[j] + delta > zero) {
					temp = b[j];
					i__1 = *n;
					for (i__ = 1; i__ <= i__1; ++i__) {
						/* L420: */
						temp -= xopt[i__] * amat[i__ + j * amat_dim1];
					}
					temp = max(temp,zero);
					if (temp >= delta) {
						temp = -temp;
					}
					rescon[j] = temp;
				}
			}
			/* L430: */
		}
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			/* L440: */
			sp[k] += sp[*npt + k];
		}

		/*     Also revise GOPT when symmetric Broyden updating is applied. */

		if (itest < 3) {
			ih = 0;
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
				i__1 = j;
				for (i__ = 1; i__ <= i__1; ++i__) {
					++ih;
					if (i__ < j) {
						gopt[j] += hq[ih] * step[i__];
					}
					/* L450: */
					gopt[i__] += hq[ih] * step[j];
				}
			}
			i__1 = *npt;
			for (k = 1; k <= i__1; ++k) {
				temp = pq[k] * sp[*npt + k];
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* L460: */
					gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
				}
			}
		}
	}

	/*     Replace the current model by the least Frobenius norm interpolant if */
	/*       this interpolant gives substantial reductions in the predictions */
	/*       of values of F at feasible points. */

	if (itest == 3) {
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			pq[k] = zero;
			/* L470: */
			w[k] = fval[k] - fval[kopt];
		}
		i__2 = nptm;
		for (j = 1; j <= i__2; ++j) {
			sum = zero;
			i__1 = *npt;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* L480: */
				sum += w[i__] * zmat[i__ + j * zmat_dim1];
			}
			if (j < idz) {
				sum = -sum;
			}
			i__1 = *npt;
			for (k = 1; k <= i__1; ++k) {
				/* L490: */
				pq[k] += sum * zmat[k + j * zmat_dim1];
			}
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			gopt[j] = zero;
			i__2 = *npt;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L500: */
				gopt[j] += w[i__] * bmat[i__ + j * bmat_dim1];
			}
		}
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			temp = pq[k] * sp[k];
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* L510: */
				gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
			}
		}
		i__1 = nh;
		for (ih = 1; ih <= i__1; ++ih) {
			/* L520: */
			hq[ih] = zero;
		}
	}

	/*     If a trust region step has provided a sufficient decrease in F, then */
	/*       branch for another trust region calculation. Every iteration that */
	/*       takes a model step is followed by an attempt to take a trust region */
	/*       step. */

	knew = 0;
	if (ksave > 0) {
		goto L20;
	}
	if (ratio >= tenth) {
		goto L20;
	}

	/*     Alternatively, find out if the interpolation points are close enough */
	/*       to the best point so far. */

L530:
	/* Computing MAX */
	d__1 = delta * delta, d__2 = rho * 4. * rho;
	distsq = max(d__1,d__2);
	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		sum = zero;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			/* L540: */
			/* Computing 2nd power */
			d__1 = xpt[k + j * xpt_dim1] - xopt[j];
			sum += d__1 * d__1;
		}
		if (sum > distsq) {
			knew = k;
			distsq = sum;
		}
		/* L550: */
	}

	/*     If KNEW is positive, then branch back for the next iteration, which */
	/*       will generate a "model step". Otherwise, if the current iteration */
	/*       has reduced F, or if DELTA was above its lower bound when the last */
	/*       trust region step was calculated, then try a "trust region" step */
	/*       instead. */

	if (knew > 0) {
		goto L20;
	}
	knew = 0;
	if (fopt < fsave) {
		goto L20;
	}
	if (delsav > rho) {
		goto L20;
	}

	/*     The calculations with the current value of RHO are complete. */
	/*       Pick the next value of RHO. */

L560:
	if (rho > *rhoend) {
		delta = half * rho;
		if (rho > *rhoend * 250.) {
			rho = tenth * rho;
		} else if (rho <= *rhoend * 16.) {
			rho = *rhoend;
		} else {
			rho = sqrt(rho * *rhoend);
		}
		delta = max(delta,rho);
		if (*iprint >= 2) {
			if (*iprint >= 3) {
				GAI::LogError("(5x)");

			}
			GAI::LogError("New RHO =\002,1pd11.4,5x,\002Number of\002,\002 function values =\002,i6");
			do_fio(&c__1, (Char *)&rho, (ftnlen)sizeof(Float64));
			do_fio(&c__1, (Char *)&nf, (ftnlen)sizeof(Int32));

			GAI::LogError("Least value of F =\002,1pd23.15,9x,\002The corresponding X is:\002/(2x,5d15.6)");
			do_fio(&c__1, (Char *)&fopt, (ftnlen)sizeof(Float64));
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				d__1 = xbase[i__] + xopt[i__];
				do_fio(&c__1, (Char *)&d__1, (ftnlen)sizeof(Float64));
			}

		}
		goto L10;
	}

	/*     Return from the calculation, after branching to label 220 for another */
	/*       Newton-Raphson step if it has not been tried before. */

	if (ksave == -1) {
		goto L220;
	}
L600:
	if (fopt <= f || ifeas == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L610: */
			x[i__] = xsav[i__];
		}
		f = fopt;
	}
	if (*iprint >= 1) {
		GAI::LogError("At the return from LINCOA\002,5x,\002Number of function values =\002,i6");
		do_fio(&c__1, (Char *)&nf, (ftnlen)sizeof(Int32));

		GAI::LogError("Least value of F =\002,1pd23.15,9x,\002The corresponding X is:\002/(2x,5d15.6)");
		do_fio(&c__1, (Char *)&f, (ftnlen)sizeof(Float64));
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			do_fio(&c__1, (Char *)&x[i__], (ftnlen)sizeof(Float64));
		}

	}
	w[1] = f;
	w[2] = (Float64) nf + half;
	return 0;
} /* lincob_ */


int GetAct(Int32 *n, Int32 *m, Float64 *amat, 
		   Float64 *b, Int32 *nact, Int32 *iact, Float64 *qfac, 
		   Float64 *rfac, Float64 *snorm, Float64 *resnew, Float64 *
		   resact, Float64 *g, Float64 *dw, Float64 *vlam, Float64 *w)
{
	/* System generated locals */
	Int32 amat_dim1, amat_offset, qfac_dim1, qfac_offset, i__1, i__2;
	Float64 d__1, d__2;

	/* Local variables */
	Int32 i__, j, k, l;
	Float64 dd;
	Int32 ic, jc, jw, jcp;
	Float64 one, sum, cval, tdel, ctol, temp, sval, cosv, zero, test, sinv, tiny;
	Int32 idiag, jdiag, iflag;
	Float64 rdiag, ddsav;
	Int32 nactp;
	Float64 dnorm, sprod, vmult, violmx;


	/*     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms */
	/*       with these names in SUBROUTINE LINCOB. The current values must be */
	/*       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when */
	/*       GETACT changes the current active set. */
	/*     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these */
	/*       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are */
	/*       also kept up to date. */
	/*     VLAM and W are used for working space, the vector VLAM being reserved */
	/*       for the Lagrange multipliers of the calculation. Their lengths must */
	/*       be at least N. */
	/*     The main purpose of GETACT is to pick the current active set. It is */
	/*       defined by the property that the projection of -G into the space */
	/*       orthogonal to the active constraint normals is as large as possible, */
	/*       subject to this projected steepest descent direction moving no closer */
	/*       to the boundary of every constraint whose current residual is at most */
	/*       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are */
	/*       all appropriate to this choice of active set. */
	/*     Occasionally this projected direction is zero, and then the final value */
	/*       of W(1) is set to zero. Otherwise, the direction itself is returned */
	/*       in DW, and W(1) is set to the square of the length of the direction. */

	/*     Set some constants and a temporary VLAM. */

	/* Parameter adjustments */
	qfac_dim1 = *n;
	qfac_offset = 1 + qfac_dim1;
	qfac -= qfac_offset;
	amat_dim1 = *n;
	amat_offset = 1 + amat_dim1;
	amat -= amat_offset;
	--b;
	--iact;
	--rfac;
	--resnew;
	--resact;
	--g;
	--dw;
	--vlam;
	--w;

	/* Function Body */
	one = 1.;
	tiny = 1e-60;
	zero = 0.;
	tdel = *snorm * .2;
	ddsav = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing 2nd power */
		d__1 = g[i__];
		ddsav += d__1 * d__1;
		/* L10: */
		vlam[i__] = zero;
	}
	ddsav += ddsav;

	/*     Set the initial QFAC to the identity matrix in the case NACT=0. */

	if (*nact == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
				/* L20: */
				qfac[i__ + j * qfac_dim1] = zero;
			}
			/* L30: */
			qfac[i__ + i__ * qfac_dim1] = one;
		}
		goto L100;
	}

	/*     Remove any constraints from the initial active set whose residuals */
	/*       exceed TDEL. */

	iflag = 1;
	ic = *nact;
L40:
	if (resact[ic] > tdel) {
		goto L800;
	}
L50:
	--ic;
	if (ic > 0) {
		goto L40;
	}

	/*     Remove any constraints from the initial active set whose Lagrange */
	/*       multipliers are nonnegative, and set the surviving multipliers. */

	iflag = 2;
L60:
	if (*nact == 0) {
		goto L100;
	}
	ic = *nact;
L70:
	temp = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* L80: */
		temp += qfac[i__ + ic * qfac_dim1] * g[i__];
	}
	idiag = (ic * ic + ic) / 2;
	if (ic < *nact) {
		jw = idiag + ic;
		i__1 = *nact;
		for (j = ic + 1; j <= i__1; ++j) {
			temp -= rfac[jw] * vlam[j];
			/* L90: */
			jw += j;
		}
	}
	if (temp >= zero) {
		goto L800;
	}
	vlam[ic] = temp / rfac[idiag];
	--ic;
	if (ic > 0) {
		goto L70;
	}

	/*     Set the new search direction D. Terminate if the 2-norm of D is zero */
	/*       or does not decrease, or if NACT=N holds. The situation NACT=N */
	/*       occurs for sufficiently large SNORM if the origin is in the convex */
	/*       hull of the constraint gradients. */

L100:
	if (*nact == *n) {
		goto L290;
	}
	i__1 = *n;
	for (j = *nact + 1; j <= i__1; ++j) {
		w[j] = zero;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L110: */
			w[j] += qfac[i__ + j * qfac_dim1] * g[i__];
		}
	}
	dd = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
		dw[i__] = zero;
		i__1 = *n;
		for (j = *nact + 1; j <= i__1; ++j) {
			/* L120: */
			dw[i__] -= w[j] * qfac[i__ + j * qfac_dim1];
		}
		/* L130: */
		/* Computing 2nd power */
		d__1 = dw[i__];
		dd += d__1 * d__1;
	}
	if (dd >= ddsav) {
		goto L290;
	}
	if (dd == zero) {
		goto L300;
	}
	ddsav = dd;
	dnorm = sqrt(dd);

	/*     Pick the next Int32 L or terminate, a positive value of L being */
	/*       the index of the most violated constraint. The purpose of CTOL */
	/*       below is to estimate whether a positive value of VIOLMX may be */
	/*       due to computer rounding errors. */

	l = 0;
	if (*m > 0) {
		test = dnorm / *snorm;
		violmx = zero;
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
			if (resnew[j] > zero && resnew[j] <= tdel) {
				sum = zero;
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L140: */
					sum += amat[i__ + j * amat_dim1] * dw[i__];
				}
				if (sum > test * resnew[j]) {
					if (sum > violmx) {
						l = j;
						violmx = sum;
					}
				}
			}
			/* L150: */
		}
		ctol = zero;
		temp = dnorm * .01;
		if (violmx > zero && violmx < temp) {
			if (*nact > 0) {
				i__2 = *nact;
				for (k = 1; k <= i__2; ++k) {
					j = iact[k];
					sum = zero;
					i__1 = *n;
					for (i__ = 1; i__ <= i__1; ++i__) {
						/* L160: */
						sum += dw[i__] * amat[i__ + j * amat_dim1];
					}
					/* L170: */
					/* Computing MAX */
					d__1 = ctol, d__2 = abs(sum);
					ctol = max(d__1,d__2);
				}
			}
		}
	}
	w[1] = one;
	if (l == 0) {
		goto L300;
	}
	if (violmx <= ctol * 10.) {
		goto L300;
	}

	/*     Apply Givens rotations to the last (N-NACT) columns of QFAC so that */
	/*       the first (NACT+1) columns of QFAC are the ones required for the */
	/*       addition of the L-th constraint, and add the appropriate column */
	/*       to RFAC. */

	nactp = *nact + 1;
	idiag = (nactp * nactp - nactp) / 2;
	rdiag = zero;
	for (j = *n; j >= 1; --j) {
		sprod = zero;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L180: */
			sprod += qfac[i__ + j * qfac_dim1] * amat[i__ + l * amat_dim1];
		}
		if (j <= *nact) {
			rfac[idiag + j] = sprod;
		} else {
			if (abs(rdiag) <= abs(sprod) * 1e-20) {
				rdiag = sprod;
			} else {
				temp = sqrt(sprod * sprod + rdiag * rdiag);
				cosv = sprod / temp;
				sinv = rdiag / temp;
				rdiag = temp;
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp = cosv * qfac[i__ + j * qfac_dim1] + sinv * qfac[i__ 
						+ (j + 1) * qfac_dim1];
					qfac[i__ + (j + 1) * qfac_dim1] = -sinv * qfac[i__ + j * 
						qfac_dim1] + cosv * qfac[i__ + (j + 1) * 
						qfac_dim1];
					/* L190: */
					qfac[i__ + j * qfac_dim1] = temp;
				}
			}
		}
		/* L200: */
	}
	if (rdiag < zero) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L210: */
			qfac[i__ + nactp * qfac_dim1] = -qfac[i__ + nactp * qfac_dim1];
		}
	}
	rfac[idiag + nactp] = abs(rdiag);
	*nact = nactp;
	iact[*nact] = l;
	resact[*nact] = resnew[l];
	vlam[*nact] = zero;
	resnew[l] = zero;

	/*     Set the components of the vector VMU in W. */

L220:
	/* Computing 2nd power */
	d__1 = rfac[(*nact * *nact + *nact) / 2];
	w[*nact] = one / (d__1 * d__1);
	if (*nact > 1) {
		for (i__ = *nact - 1; i__ >= 1; --i__) {
			idiag = (i__ * i__ + i__) / 2;
			jw = idiag + i__;
			sum = zero;
			i__2 = *nact;
			for (j = i__ + 1; j <= i__2; ++j) {
				sum -= rfac[jw] * w[j];
				/* L230: */
				jw += j;
			}
			/* L240: */
			w[i__] = sum / rfac[idiag];
		}
	}

	/*     Calculate the multiple of VMU to subtract from VLAM, and update VLAM. */

	vmult = violmx;
	ic = 0;
	j = 1;
L250:
	if (j < *nact) {
		if (vlam[j] >= vmult * w[j]) {
			ic = j;
			vmult = vlam[j] / w[j];
		}
		++j;
		goto L250;
	}
	i__2 = *nact;
	for (j = 1; j <= i__2; ++j) {
		/* L260: */
		vlam[j] -= vmult * w[j];
	}
	if (ic > 0) {
		vlam[ic] = zero;
	}
	/* Computing MAX */
	d__1 = violmx - vmult;
	violmx = max(d__1,zero);
	if (ic == 0) {
		violmx = zero;
	}

	/*     Reduce the active set if necessary, so that all components of the */
	/*       new VLAM are negative, with resetting of the residuals of the */
	/*       constraints that become inactive. */

	iflag = 3;
	ic = *nact;
L270:
	if (vlam[ic] < zero) {
		goto L280;
	}
	/* Computing MAX */
	d__1 = resact[ic];
	resnew[iact[ic]] = max(d__1,tiny);
	goto L800;
L280:
	--ic;
	if (ic > 0) {
		goto L270;
	}

	/*     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds, */
	/*       as then the active constraints imply D=0. Otherwise, go to label */
	/*       100, to calculate the new D and to test for termination. */

	if (violmx > zero) {
		goto L220;
	}
	if (*nact < *n) {
		goto L100;
	}
L290:
	dd = zero;
L300:
	w[1] = dd;
	return 0;

	/*     These instructions rearrange the active constraints so that the new */
	/*       value of IACT(NACT) is the old value of IACT(IC). A sequence of */
	/*       Givens rotations is applied to the current QFAC and RFAC. Then NACT */
	/*       is reduced by one. */

L800:
	/* Computing MAX */
	d__1 = resact[ic];
	resnew[iact[ic]] = max(d__1,tiny);
	jc = ic;
L810:
	if (jc < *nact) {
		jcp = jc + 1;
		idiag = jc * jcp / 2;
		jw = idiag + jcp;
		/* Computing 2nd power */
		d__1 = rfac[jw - 1];
		/* Computing 2nd power */
		d__2 = rfac[jw];
		temp = sqrt(d__1 * d__1 + d__2 * d__2);
		cval = rfac[jw] / temp;
		sval = rfac[jw - 1] / temp;
		rfac[jw - 1] = sval * rfac[idiag];
		rfac[jw] = cval * rfac[idiag];
		rfac[idiag] = temp;
		if (jcp < *nact) {
			i__2 = *nact;
			for (j = jcp + 1; j <= i__2; ++j) {
				temp = sval * rfac[jw + jc] + cval * rfac[jw + jcp];
				rfac[jw + jcp] = cval * rfac[jw + jc] - sval * rfac[jw + jcp];
				rfac[jw + jc] = temp;
				/* L820: */
				jw += j;
			}
		}
		jdiag = idiag - jc;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			if (i__ < jc) {
				temp = rfac[idiag + i__];
				rfac[idiag + i__] = rfac[jdiag + i__];
				rfac[jdiag + i__] = temp;
			}
			temp = sval * qfac[i__ + jc * qfac_dim1] + cval * qfac[i__ + jcp *
				qfac_dim1];
			qfac[i__ + jcp * qfac_dim1] = cval * qfac[i__ + jc * qfac_dim1] - 
				sval * qfac[i__ + jcp * qfac_dim1];
			/* L830: */
			qfac[i__ + jc * qfac_dim1] = temp;
		}
		iact[jc] = iact[jcp];
		resact[jc] = resact[jcp];
		vlam[jc] = vlam[jcp];
		jc = jcp;
		goto L810;
	}
	--(*nact);
	switch (iflag) {
	case 1:  goto L50;
	case 2:  goto L60;
	case 3:  goto L280;
	}
	return 0;
} /* getact_ */



template<typename tFunc> int Linco_Initialise(tFunc &lpFunction,Int32 *n, Int32 *npt, Int32 *m, Float64 
											  *amat, Float64 *b, Float64 *x, Float64 *rhobeg, Int32 *
											  iprint, Float64 *xbase, Float64 *xpt, Float64 *fval, 
											  Float64 *xsav, Float64 *xopt, Float64 *gopt, Int32 *kopt, 
											  Float64 *hq, Float64 *pq, Float64 *bmat, Float64 *zmat, 
											  Int32 *idz, Int32 *ndim, Float64 *sp, Float64 *rescon, 
											  Float64 *step, Float64 *pqw, Float64 *w)
{
	
	/* System generated locals */
	Int32 amat_dim1, amat_offset, xpt_dim1, xpt_offset, bmat_dim1, 
		bmat_offset, zmat_dim1, zmat_offset, i__1, i__2, i__3;


	/* Local variables */
	Float64 f;
	Int32 i__, j, k, nf, jp;
	Float64 one;
	Int32 ipt, jpt;
	Float64 half, feas, bigv;
	Int32 jsav;
	Float64 temp;
	Int32 nptm;
	Float64 zero, test;
	Int32 kbase;
	Float64 recip, reciq, resid;
	Int32 itemp;
	Float64 rhosq;


	/*     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL, */
	/*       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the */
	/*       same as the corresponding arguments in SUBROUTINE LINCOB. */
	/*     KOPT is set to the Int32 such that XPT(KOPT,.) is the initial trust */
	/*       region centre. */
	/*     IDZ is going to be set to one, so that every element of Diag(DZ) is */
	/*       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the */
	/*       factorization of the leading NPT by NPT submatrix of H. */
	/*     STEP, PQW and W are used for working space, the arrays STEP and PQW */
	/*       being taken from LINCOB. The length of W must be at least N+NPT. */

	/*     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT */
	/*       for the first iteration, an important feature being that, if any of */
	/*       of the columns of XPT is an infeasible point, then the largest of */
	/*       the constraint violations there is at least 0.2*RHOBEG. It also sets */
	/*       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON. */

	/*     Set some constants. */

	/* Parameter adjustments */
	amat_dim1 = *n;
	amat_offset = 1 + amat_dim1;
	amat -= amat_offset;
	zmat_dim1 = *npt;
	zmat_offset = 1 + zmat_dim1;
	zmat -= zmat_offset;
	xpt_dim1 = *npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	--b;
	--x;
	--xbase;
	--fval;
	--xsav;
	--xopt;
	--gopt;
	--hq;
	--pq;
	bmat_dim1 = *ndim;
	bmat_offset = 1 + bmat_dim1;
	bmat -= bmat_offset;
	--sp;
	--rescon;
	--step;
	--pqw;
	--w;

	/* Function Body */
	half = .5;
	one = 1.;
	zero = 0.;
	nptm = *npt - *n - 1;
	rhosq = *rhobeg * *rhobeg;
	recip = one / rhosq;
	reciq = sqrt(half) / rhosq;
	test = *rhobeg * .2;
	*idz = 1;
	kbase = 1;

	/*     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		xbase[j] = x[j];
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			/* L10: */
			xpt[k + j * xpt_dim1] = zero;
		}
		i__2 = *ndim;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L20: */
			bmat[i__ + j * bmat_dim1] = zero;
		}
	}
	i__2 = *npt;
	for (k = 1; k <= i__2; ++k) {
		sp[k] = zero;
		i__1 = *npt - *n - 1;
		for (j = 1; j <= i__1; ++j) {
			/* L30: */
			zmat[k + j * zmat_dim1] = zero;
		}
	}

	/*     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT], */
	/*       but they may be altered later to make a constraint violation */
	/*       sufficiently large. The initial nonzero elements of BMAT and of */
	/*       the first min[N,NPT-N-1] columns of ZMAT are set also. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		xpt[j + 1 + j * xpt_dim1] = *rhobeg;
		if (j < *npt - *n) {
			jp = *n + j + 1;
			xpt[jp + j * xpt_dim1] = -(*rhobeg);
			bmat[j + 1 + j * bmat_dim1] = half / *rhobeg;
			bmat[jp + j * bmat_dim1] = -half / *rhobeg;
			zmat[j * zmat_dim1 + 1] = -reciq - reciq;
			zmat[j + 1 + j * zmat_dim1] = reciq;
			zmat[jp + j * zmat_dim1] = reciq;
		} else {
			bmat[j * bmat_dim1 + 1] = -one / *rhobeg;
			bmat[j + 1 + j * bmat_dim1] = one / *rhobeg;
			bmat[*npt + j + j * bmat_dim1] = -half * rhosq;
		}
		/* L40: */
	}

	/*     Set the remaining initial nonzero elements of XPT and ZMAT when the */
	/*       number of interpolation points exceeds 2*N+1. */

	if (*npt > (*n << 1) + 1) {
		i__1 = *npt - *n - 1;
		for (k = *n + 1; k <= i__1; ++k) {
			itemp = (k - 1) / *n;
			ipt = k - itemp * *n;
			jpt = ipt + itemp;
			if (jpt > *n) {
				jpt -= *n;
			}
			xpt[*n + k + 1 + ipt * xpt_dim1] = *rhobeg;
			xpt[*n + k + 1 + jpt * xpt_dim1] = *rhobeg;
			zmat[k * zmat_dim1 + 1] = recip;
			zmat[ipt + 1 + k * zmat_dim1] = -recip;
			zmat[jpt + 1 + k * zmat_dim1] = -recip;
			/* L50: */
			zmat[*n + k + 1 + k * zmat_dim1] = recip;
		}
	}

	/*     Update the constraint right hand sides to allow for the shift XBASE. */

	if (*m > 0) {
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
			temp = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L60: */
				temp += amat[i__ + j * amat_dim1] * xbase[i__];
			}
			/* L70: */
			b[j] -= temp;
		}
	}

	/*     Go through the initial points, shifting every infeasible point if */
	/*       necessary so that its constraint violation is at least 0.2*RHOBEG. */

	i__1 = *npt;
	for (nf = 1; nf <= i__1; ++nf) {
		feas = one;
		bigv = zero;
		j = 0;
L80:
		++j;
		if (j <= *m && nf >= 2) {
			resid = -b[j];
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L90: */
				resid += xpt[nf + i__ * xpt_dim1] * amat[i__ + j * amat_dim1];
			}
			if (resid <= bigv) {
				goto L80;
			}
			bigv = resid;
			jsav = j;
			if (resid <= test) {
				feas = -one;
				goto L80;
			}
			feas = zero;
		}
		if (feas < zero) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L100: */
				step[i__] = xpt[nf + i__ * xpt_dim1] + (test - bigv) * amat[
					i__ + jsav * amat_dim1];
			}
			i__2 = *npt;
			for (k = 1; k <= i__2; ++k) {
				sp[*npt + k] = zero;
				i__3 = *n;
				for (j = 1; j <= i__3; ++j) {
					/* L110: */
					sp[*npt + k] += xpt[k + j * xpt_dim1] * step[j];
				}
			}
			update_(n, npt, &xpt[xpt_offset], &bmat[bmat_offset], &zmat[
				zmat_offset], idz, ndim, &sp[1], &step[1], &kbase, &nf, &
					pqw[1], &w[1]);
				i__3 = *n;
				for (i__ = 1; i__ <= i__3; ++i__) {
					/* L120: */
					xpt[nf + i__ * xpt_dim1] = step[i__];
				}
		}

		/*     Calculate the objective function at the current interpolation point, */
		/*       and set KOPT to the index of the first trust region centre. */

		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
			/* L130: */
			x[j] = xbase[j] + xpt[nf + j * xpt_dim1];
		}
		f = feas;
		lpFunction(n, &x[1], &f);
		if (*iprint == 3) {
			s_wsfe(&io___24);
			do_fio(&c__1, (Char *)&nf, (ftnlen)sizeof(Int32));
			do_fio(&c__1, (Char *)&f, (ftnlen)sizeof(Float64));
			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
				do_fio(&c__1, (Char *)&x[i__], (ftnlen)sizeof(Float64));
			}

		}
		if (nf == 1) {
			*kopt = 1;
		} else if (f < fval[*kopt] && feas > zero) {
			*kopt = nf;
		}
		/* L150: */
		fval[nf] = f;
	}

	/*     Set PQ for the first quadratic model. */

	i__1 = nptm;
	for (j = 1; j <= i__1; ++j) {
		w[j] = zero;
		i__3 = *npt;
		for (k = 1; k <= i__3; ++k) {
			/* L160: */
			w[j] += zmat[k + j * zmat_dim1] * fval[k];
		}
	}
	i__3 = *npt;
	for (k = 1; k <= i__3; ++k) {
		pq[k] = zero;
		i__1 = nptm;
		for (j = 1; j <= i__1; ++j) {
			/* L170: */
			pq[k] += zmat[k + j * zmat_dim1] * w[j];
		}
	}

	/*     Set XOPT, SP, GOPT and HQ for the first quadratic model. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		xopt[j] = xpt[*kopt + j * xpt_dim1];
		xsav[j] = xbase[j] + xopt[j];
		/* L180: */
		gopt[j] = zero;
	}
	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		sp[k] = zero;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
			/* L190: */
			sp[k] += xpt[k + j * xpt_dim1] * xopt[j];
		}
		temp = pq[k] * sp[k];
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
			/* L200: */
			gopt[j] = gopt[j] + fval[k] * bmat[k + j * bmat_dim1] + temp * 
				xpt[k + j * xpt_dim1];
		}
	}
	i__3 = (*n * *n + *n) / 2;
	for (i__ = 1; i__ <= i__3; ++i__) {
		/* L210: */
		hq[i__] = zero;
	}

	/*     Set the initial elements of RESCON. */

	i__3 = *m;
	for (j = 1; j <= i__3; ++j) {
		temp = b[j];
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L220: */
			temp -= xopt[i__] * amat[i__ + j * amat_dim1];
		}
		temp = max(temp,zero);
		if (temp >= *rhobeg) {
			temp = -temp;
		}
		/* L230: */
		rescon[j] = temp;
	}
	return 0;
} /* prelim_ */



int trstep_(Int32 *n, Int32 *npt, Int32 *m, Float64 
							 *amat, Float64 *b, Float64 *xpt, Float64 *hq, Float64 *pq,
							 Int32 *nact, Int32 *iact, Float64 *rescon, Float64 *qfac, 
							 Float64 *rfac, Float64 *snorm, Float64 *step, Float64 *g, 
							 Float64 *resnew, Float64 *resact, Float64 *d__, Float64 *
							 dw, Float64 *w)
{
	/* System generated locals */
	Int32 amat_dim1, amat_offset, xpt_dim1, xpt_offset, qfac_dim1, 
		qfac_offset, i__1, i__2;
	Float64 d__1, d__2;

	/* Local variables */
	Int32 i__, j, k;
	Float64 ad, dd, dg;
	Int32 ih;
	Float64 ds;
	Int32 ir;
	Float64 ss, dgd, adw, one, wgd, rhs, sum, half, beta;
	Int32 jsav;
	Float64 temp, zero, tiny, snsq, gamma, alpbd, alpha, scale;
	Int32 ncall;
	Float64 alphm, alpht, ctest;
	Float64 reduct, resmax;
	Int32 icount;


	/*     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC */
	/*       are the same as the terms with these names in LINCOB. If RESCON(J) */
	/*       is negative, then |RESCON(J)| must be no less than the trust region */
	/*       radius, so that the J-th constraint can be ignored. */
	/*     SNORM is set to the trust region radius DELTA initially. On the */
	/*       return, however, it is the length of the calculated STEP, which is */
	/*       set to zero if the constraints do not allow a long enough step. */
	/*     STEP is the total calculated step so far from the trust region centre, */
	/*       its final value being given by the sequence of CG iterations, which */
	/*       terminate if the trust region boundary is reached. */
	/*     G must be set on entry to the gradient of the quadratic model at the */
	/*       trust region centre. It is used as working space, however, and is */
	/*       always the gradient of the model at the current STEP, except that */
	/*       on return the value of G(1) is set to ONE instead of to ZERO if */
	/*       and only if GETACT is called more than once. */
	/*     RESNEW, RESACT, D, DW and W are used for working space. A negative */
	/*       value of RESNEW(J) indicates that the J-th constraint does not */
	/*       restrict the CG steps of the current trust region calculation, a */
	/*       zero value of RESNEW(J) indicates that the J-th constraint is active, */
	/*       and otherwise RESNEW(J) is set to the greater of TINY and the actual */
	/*       residual of the J-th constraint for the current STEP. RESACT holds */
	/*       the residuals of the active constraints, which may be positive. */
	/*       D is the search direction of each line search. DW is either another */
	/*       search direction or the change in gradient along D. The length of W */
	/*       must be at least MAX[M,2*N]. */

	/*     Set some numbers for the conjugate gradient iterations. */

	/* Parameter adjustments */
	qfac_dim1 = *n;
	qfac_offset = 1 + qfac_dim1;
	qfac -= qfac_offset;
	amat_dim1 = *n;
	amat_offset = 1 + amat_dim1;
	amat -= amat_offset;
	xpt_dim1 = *npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	--b;
	--hq;
	--pq;
	--iact;
	--rescon;
	--rfac;
	--step;
	--g;
	--resnew;
	--resact;
	--d__;
	--dw;
	--w;

	/* Function Body */
	half = .5;
	one = 1.;
	tiny = 1e-60;
	zero = 0.;
	ctest = .01;
	snsq = *snorm * *snorm;

	/*     Set the initial elements of RESNEW, RESACT and STEP. */

	if (*m > 0) {
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
			resnew[j] = rescon[j];
			if (rescon[j] >= *snorm) {
				resnew[j] = -one;
			} else if (rescon[j] >= zero) {
				/* Computing MAX */
				d__1 = resnew[j];
				resnew[j] = max(d__1,tiny);
			}
			/* L10: */
		}
		if (*nact > 0) {
			i__1 = *nact;
			for (k = 1; k <= i__1; ++k) {
				resact[k] = rescon[iact[k]];
				/* L20: */
				resnew[iact[k]] = zero;
			}
		}
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* L30: */
		step[i__] = zero;
	}
	ss = zero;
	reduct = zero;
	ncall = 0;

	/*     GETACT picks the active set for the current STEP. It also sets DW to */
	/*       the vector closest to -G that is orthogonal to the normals of the */
	/*       active constraints. DW is scaled to have length 0.2*SNORM, as then */
	/*       a move of DW from STEP is allowed by the linear constraints. */

L40:
	++ncall;
	GetAct(n, m, &amat[amat_offset], &b[1], nact, &iact[1], &qfac[
		qfac_offset], &rfac[1], snorm, &resnew[1], &resact[1], &g[1], &dw[
			1], &w[1], &w[*n + 1]);
			if (w[*n + 1] == zero) {
				goto L320;
			}
			scale = *snorm * .2 / sqrt(w[*n + 1]);
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* L50: */
				dw[i__] = scale * dw[i__];
			}

			/*     If the modulus of the residual of an active constraint is substantial, */
			/*       then set D to the shortest move from STEP to the boundaries of the */
			/*       active constraints. */

			resmax = zero;
			if (*nact > 0) {
				i__1 = *nact;
				for (k = 1; k <= i__1; ++k) {
					/* L60: */
					/* Computing MAX */
					d__1 = resmax, d__2 = resact[k];
					resmax = max(d__1,d__2);
				}
			}
			gamma = zero;
			if (resmax > *snorm * 1e-4) {
				ir = 0;
				i__1 = *nact;
				for (k = 1; k <= i__1; ++k) {
					temp = resact[k];
					if (k >= 2) {
						i__2 = k - 1;
						for (i__ = 1; i__ <= i__2; ++i__) {
							++ir;
							/* L70: */
							temp -= rfac[ir] * w[i__];
						}
					}
					++ir;
					/* L80: */
					w[k] = temp / rfac[ir];
				}
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					d__[i__] = zero;
					i__2 = *nact;
					for (k = 1; k <= i__2; ++k) {
						/* L90: */
						d__[i__] += w[k] * qfac[i__ + k * qfac_dim1];
					}
				}

				/*     The vector D that has just been calculated is also the shortest move */
				/*       from STEP+DW to the boundaries of the active constraints. Set GAMMA */
				/*       to the greatest steplength of this move that satisfies the trust */
				/*       region bound. */

				rhs = snsq;
				ds = zero;
				dd = zero;
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					sum = step[i__] + dw[i__];
					rhs -= sum * sum;
					ds += d__[i__] * sum;
					/* L100: */
					/* Computing 2nd power */
					d__1 = d__[i__];
					dd += d__1 * d__1;
				}
				if (rhs > zero) {
					temp = sqrt(ds * ds + dd * rhs);
					if (ds <= zero) {
						gamma = (temp - ds) / dd;
					} else {
						gamma = rhs / (temp + ds);
					}
				}

				/*     Reduce the steplength GAMMA if necessary so that the move along D */
				/*       also satisfies the linear constraints. */

				j = 0;
L110:
				if (gamma > zero) {
					++j;
					if (resnew[j] > zero) {
						ad = zero;
						adw = zero;
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							ad += amat[i__ + j * amat_dim1] * d__[i__];
							/* L120: */
							adw += amat[i__ + j * amat_dim1] * dw[i__];
						}
						if (ad > zero) {
							/* Computing MAX */
							d__1 = (resnew[j] - adw) / ad;
							temp = max(d__1,zero);
							gamma = min(gamma,temp);
						}
					}
					if (j < *m) {
						goto L110;
					}
				}
				gamma = min(gamma,one);
			}

			/*     Set the next direction for seeking a reduction in the model function */
			/*       subject to the trust region bound and the linear constraints. */

			if (gamma <= zero) {
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* L130: */
					d__[i__] = dw[i__];
				}
				icount = *nact;
			} else {
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* L140: */
					d__[i__] = dw[i__] + gamma * d__[i__];
				}
				icount = *nact - 1;
			}
			alpbd = one;

			/*     Set ALPHA to the steplength from STEP along D to the trust region */
			/*       boundary. Return if the first derivative term of this step is */
			/*       sufficiently small or if no further progress is possible. */

L150:
			++icount;
			rhs = snsq - ss;
			if (rhs <= zero) {
				goto L320;
			}
			dg = zero;
			ds = zero;
			dd = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				dg += d__[i__] * g[i__];
				ds += d__[i__] * step[i__];
				/* L160: */
				/* Computing 2nd power */
				d__1 = d__[i__];
				dd += d__1 * d__1;
			}
			if (dg >= zero) {
				goto L320;
			}
			temp = sqrt(rhs * dd + ds * ds);
			if (ds <= zero) {
				alpha = (temp - ds) / dd;
			} else {
				alpha = rhs / (temp + ds);
			}
			if (-alpha * dg <= ctest * reduct) {
				goto L320;
			}

			/*     Set DW to the change in gradient along D. */

			ih = 0;
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
				dw[j] = zero;
				i__1 = j;
				for (i__ = 1; i__ <= i__1; ++i__) {
					++ih;
					if (i__ < j) {
						dw[j] += hq[ih] * d__[i__];
					}
					/* L170: */
					dw[i__] += hq[ih] * d__[j];
				}
			}
			i__1 = *npt;
			for (k = 1; k <= i__1; ++k) {
				temp = zero;
				i__2 = *n;
				for (j = 1; j <= i__2; ++j) {
					/* L180: */
					temp += xpt[k + j * xpt_dim1] * d__[j];
				}
				temp = pq[k] * temp;
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* L190: */
					dw[i__] += temp * xpt[k + i__ * xpt_dim1];
				}
			}

			/*     Set DGD to the curvature of the model along D. Then reduce ALPHA if */
			/*       necessary to the value that minimizes the model. */

			dgd = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* L200: */
				dgd += d__[i__] * dw[i__];
			}
			alpht = alpha;
			if (dg + alpha * dgd > zero) {
				alpha = -dg / dgd;
			}

			/*     Make a further reduction in ALPHA if necessary to preserve feasibility, */
			/*       and put some scalar products of D with constraint gradients in W. */

			alphm = alpha;
			jsav = 0;
			if (*m > 0) {
				i__2 = *m;
				for (j = 1; j <= i__2; ++j) {
					ad = zero;
					if (resnew[j] > zero) {
						i__1 = *n;
						for (i__ = 1; i__ <= i__1; ++i__) {
							/* L210: */
							ad += amat[i__ + j * amat_dim1] * d__[i__];
						}
						if (alpha * ad > resnew[j]) {
							alpha = resnew[j] / ad;
							jsav = j;
						}
					}
					/* L220: */
					w[j] = ad;
				}
			}
			alpha = max(alpha,alpbd);
			alpha = min(alpha,alphm);
			if (icount == *nact) {
				alpha = min(alpha,one);
			}

			/*     Update STEP, G, RESNEW, RESACT and REDUCT. */

			ss = zero;
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
				step[i__] += alpha * d__[i__];
				/* Computing 2nd power */
				d__1 = step[i__];
				ss += d__1 * d__1;
				/* L230: */
				g[i__] += alpha * dw[i__];
			}
			if (*m > 0) {
				i__2 = *m;
				for (j = 1; j <= i__2; ++j) {
					if (resnew[j] > zero) {
						/* Computing MAX */
						d__1 = resnew[j] - alpha * w[j];
						resnew[j] = max(d__1,tiny);
					}
					/* L240: */
				}
			}
			if (icount == *nact && *nact > 0) {
				i__2 = *nact;
				for (k = 1; k <= i__2; ++k) {
					/* L250: */
					resact[k] = (one - gamma) * resact[k];
				}
			}
			reduct -= alpha * (dg + half * alpha * dgd);

			/*     Test for termination. Branch to label 40 if there is a new active */
			/*       constraint and if the distance from STEP to the trust region */
			/*       boundary is at least 0.2*SNORM. */

			if (alpha == alpht) {
				goto L320;
			}
			temp = -alphm * (dg + half * alphm * dgd);
			if (temp <= ctest * reduct) {
				goto L320;
			}
			if (jsav > 0) {
				if (ss <= snsq * .64) {
					goto L40;
				}
				goto L320;
			}
			if (icount == *n) {
				goto L320;
			}

			/*     Calculate the next search direction, which is conjugate to the */
			/*       previous one except in the case ICOUNT=NACT. */

			if (*nact > 0) {
				i__2 = *n;
				for (j = *nact + 1; j <= i__2; ++j) {
					w[j] = zero;
					i__1 = *n;
					for (i__ = 1; i__ <= i__1; ++i__) {
						/* L260: */
						w[j] += g[i__] * qfac[i__ + j * qfac_dim1];
					}
				}
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					temp = zero;
					i__2 = *n;
					for (j = *nact + 1; j <= i__2; ++j) {
						/* L270: */
						temp += qfac[i__ + j * qfac_dim1] * w[j];
					}
					/* L280: */
					w[*n + i__] = temp;
				}
			} else {
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L290: */
					w[*n + i__] = g[i__];
				}
			}
			if (icount == *nact) {
				beta = zero;
			} else {
				wgd = zero;
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L300: */
					wgd += w[*n + i__] * dw[i__];
				}
				beta = wgd / dgd;
			}
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* L310: */
				d__[i__] = -w[*n + i__] + beta * d__[i__];
			}
			alpbd = zero;
			goto L150;

			/*     Return from the subroutine. */

L320:
			*snorm = zero;
			if (reduct > zero) {
				*snorm = sqrt(ss);
			}
			g[1] = zero;
			if (ncall > 1) {
				g[1] = one;
			}
			return 0;
} /* trstep_ */


/* Subroutine */ int update_(Int32 *n, Int32 *npt, Float64 *xpt, 
							 Float64 *bmat, Float64 *zmat, Int32 *idz, Int32 *ndim, 
							 Float64 *sp, Float64 *step, Int32 *kopt, Int32 *knew, 
							 Float64 *vlag, Float64 *w)
{
	/* System generated locals */
	Int32 xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1, 
		zmat_offset, i__1, i__2;
	Float64 d__1, d__2;

	
	/* Local variables */
	Int32 i__, j, k, ja, jb, jl, jp;
	Float64 dx, one, tau, sum, ssq, half, beta, temp, bsum;
	Int32 nptm;
	Float64 zero, hdiag;
	Int32 iflag;
	Float64 scala, scalb, alpha, denom, tempa, tempb, tausq, denabs,
		denmax, distsq, sqrtdn;


	/*     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are */
	/*       identical to the corresponding arguments in SUBROUTINE LINCOB. */
	/*     KOPT is such that XPT(KOPT,.) is the current trust region centre. */
	/*     KNEW on exit is usually positive, and then it is the index of an */
	/*       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.). */
	/*       It is set on entry either to its final value or to 0. In the latter */
	/*       case, the final value of KNEW is chosen to maximize the denominator */
	/*       of the matrix updating formula times a weighting factor. */
	/*     VLAG and W are used for working space, the first NPT+N elements of */
	/*       both of these vectors being required. */

	/*     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being */
	/*       the ones that are suitable after the shift of the KNEW-th point to */
	/*       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero */
	/*       occurs if the calculation fails due to a zero denominator in the */
	/*       updating formula, which should never happen. */

	/*     Set some constants. */

	/* Parameter adjustments */
	zmat_dim1 = *npt;
	zmat_offset = 1 + zmat_dim1;
	zmat -= zmat_offset;
	xpt_dim1 = *npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	bmat_dim1 = *ndim;
	bmat_offset = 1 + bmat_dim1;
	bmat -= bmat_offset;
	--sp;
	--step;
	--vlag;
	--w;

	/* Function Body */
	half = .5;
	one = 1.;
	zero = 0.;
	nptm = *npt - *n - 1;

	/*     Calculate VLAG and BETA for the current choice of STEP. The first NPT */
	/*       elements of VLAG are set to the values of the Lagrange functions at */
	/*       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held */
	/*       in W, where W_check is defined in a paper on the updating method. */

	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		w[k] = sp[*npt + k] * (half * sp[*npt + k] + sp[k]);
		sum = zero;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			/* L10: */
			sum += bmat[k + j * bmat_dim1] * step[j];
		}
		/* L20: */
		vlag[k] = sum;
	}
	beta = zero;
	i__1 = nptm;
	for (k = 1; k <= i__1; ++k) {
		sum = zero;
		i__2 = *npt;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L30: */
			sum += zmat[i__ + k * zmat_dim1] * w[i__];
		}
		if (k < *idz) {
			beta += sum * sum;
			sum = -sum;
		} else {
			beta -= sum * sum;
		}
		i__2 = *npt;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L40: */
			vlag[i__] += sum * zmat[i__ + k * zmat_dim1];
		}
	}
	bsum = zero;
	dx = zero;
	ssq = zero;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
		sum = zero;
		i__1 = *npt;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L50: */
			sum += w[i__] * bmat[i__ + j * bmat_dim1];
		}
		bsum += sum * step[j];
		jp = *npt + j;
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
			/* L60: */
			sum += bmat[jp + k * bmat_dim1] * step[k];
		}
		vlag[jp] = sum;
		bsum += sum * step[j];
		dx += step[j] * xpt[*kopt + j * xpt_dim1];
		/* L70: */
		/* Computing 2nd power */
		d__1 = step[j];
		ssq += d__1 * d__1;
	}
	beta = dx * dx + ssq * (sp[*kopt] + dx + dx + half * ssq) + beta - bsum;
	vlag[*kopt] += one;

	/*     If KNEW is zero initially, then pick the index of the interpolation */
	/*       point to be deleted, by maximizing the absolute value of the */
	/*       denominator of the updating formula times a weighting factor. */


	if (*knew == 0) {
		denmax = zero;
		i__2 = *npt;
		for (k = 1; k <= i__2; ++k) {
			hdiag = zero;
			i__1 = nptm;
			for (j = 1; j <= i__1; ++j) {
				temp = one;
				if (j < *idz) {
					temp = -one;
				}
				/* L80: */
				/* Computing 2nd power */
				d__1 = zmat[k + j * zmat_dim1];
				hdiag += temp * (d__1 * d__1);
			}
			/* Computing 2nd power */
			d__2 = vlag[k];
			denabs = (d__1 = beta * hdiag + d__2 * d__2, abs(d__1));
			distsq = zero;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				/* L90: */
				/* Computing 2nd power */
				d__1 = xpt[k + j * xpt_dim1] - xpt[*kopt + j * xpt_dim1];
				distsq += d__1 * d__1;
			}
			temp = denabs * distsq * distsq;
			if (temp > denmax) {
				denmax = temp;
				*knew = k;
			}
			/* L100: */
		}
	}

	/*     Apply the rotations that put zeros in the KNEW-th row of ZMAT. */

	jl = 1;
	if (nptm >= 2) {
		i__2 = nptm;
		for (j = 2; j <= i__2; ++j) {
			if (j == *idz) {
				jl = *idz;
			} else if (zmat[*knew + j * zmat_dim1] != zero) {
				/* Computing 2nd power */
				d__1 = zmat[*knew + jl * zmat_dim1];
				/* Computing 2nd power */
				d__2 = zmat[*knew + j * zmat_dim1];
				temp = sqrt(d__1 * d__1 + d__2 * d__2);
				tempa = zmat[*knew + jl * zmat_dim1] / temp;
				tempb = zmat[*knew + j * zmat_dim1] / temp;
				i__1 = *npt;
				for (i__ = 1; i__ <= i__1; ++i__) {
					temp = tempa * zmat[i__ + jl * zmat_dim1] + tempb * zmat[
						i__ + j * zmat_dim1];
						zmat[i__ + j * zmat_dim1] = tempa * zmat[i__ + j * 
							zmat_dim1] - tempb * zmat[i__ + jl * zmat_dim1];
						/* L110: */
						zmat[i__ + jl * zmat_dim1] = temp;
				}
				zmat[*knew + j * zmat_dim1] = zero;
			}
			/* L120: */
		}
	}

	/*     Put the first NPT components of the KNEW-th column of the Z Z^T matrix */
	/*       into W, and calculate the parameters of the updating formula. */

	tempa = zmat[*knew + zmat_dim1];
	if (*idz >= 2) {
		tempa = -tempa;
	}
	if (jl > 1) {
		tempb = zmat[*knew + jl * zmat_dim1];
	}
	i__2 = *npt;
	for (i__ = 1; i__ <= i__2; ++i__) {
		w[i__] = tempa * zmat[i__ + zmat_dim1];
		if (jl > 1) {
			w[i__] += tempb * zmat[i__ + jl * zmat_dim1];
		}
		/* L130: */
	}
	alpha = w[*knew];
	tau = vlag[*knew];
	tausq = tau * tau;
	denom = alpha * beta + tausq;
	vlag[*knew] -= one;
	if (denom == zero) {
		*knew = 0;
		goto L180;
	}
	sqrtdn = sqrt((abs(denom)));

	/*     Complete the updating of ZMAT when there is only one nonzero element */
	/*       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when */
	/*       the value of IDZ is going to be reduced. */

	iflag = 0;
	if (jl == 1) {
		tempa = tau / sqrtdn;
		tempb = zmat[*knew + zmat_dim1] / sqrtdn;
		i__2 = *npt;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L140: */
			zmat[i__ + zmat_dim1] = tempa * zmat[i__ + zmat_dim1] - tempb * 
				vlag[i__];
		}
		if (denom < zero) {
			if (*idz == 1) {
				*idz = 2;
			} else {
				iflag = 1;
			}
		}
	} else {

		/*     Complete the updating of ZMAT in the alternative case. */

		ja = 1;
		if (beta >= zero) {
			ja = jl;
		}
		jb = jl + 1 - ja;
		temp = zmat[*knew + jb * zmat_dim1] / denom;
		tempa = temp * beta;
		tempb = temp * tau;
		temp = zmat[*knew + ja * zmat_dim1];
		scala = one / sqrt(abs(beta) * temp * temp + tausq);
		scalb = scala * sqrtdn;
		i__2 = *npt;
		for (i__ = 1; i__ <= i__2; ++i__) {
			zmat[i__ + ja * zmat_dim1] = scala * (tau * zmat[i__ + ja * 
				zmat_dim1] - temp * vlag[i__]);
			/* L150: */
			zmat[i__ + jb * zmat_dim1] = scalb * (zmat[i__ + jb * zmat_dim1] 
			- tempa * w[i__] - tempb * vlag[i__]);
		}
		if (denom <= zero) {
			if (beta < zero) {
				++(*idz);
			} else {
				iflag = 1;
			}
		}
	}

	/*     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times */
	/*       ZMAT^T factorization gains another positive element. Then exchange */
	/*       the first and IDZ-th columns of ZMAT. */

	if (iflag == 1) {
		--(*idz);
		i__2 = *npt;
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = zmat[i__ + zmat_dim1];
			zmat[i__ + zmat_dim1] = zmat[i__ + *idz * zmat_dim1];
			/* L160: */
			zmat[i__ + *idz * zmat_dim1] = temp;
		}
	}

	/*     Finally, update the matrix BMAT. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
		jp = *npt + j;
		w[jp] = bmat[*knew + j * bmat_dim1];
		tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
		tempb = (-beta * w[jp] - tau * vlag[jp]) / denom;
		i__1 = jp;
		for (i__ = 1; i__ <= i__1; ++i__) {
			bmat[i__ + j * bmat_dim1] = bmat[i__ + j * bmat_dim1] + tempa * 
				vlag[i__] + tempb * w[i__];
			if (i__ > *npt) {
				bmat[jp + (i__ - *npt) * bmat_dim1] = bmat[i__ + j * 
					bmat_dim1];
			}
			/* L170: */
		}
	}
L180:
	return 0;
} /* update_ */



int qmstep_(Int32 *n, Int32 *npt, Int32 *m, Float64 
							 *amat, Float64 *b, Float64 *xpt, Float64 *xopt, Int32 *
							 nact, Int32 *iact, Float64 *rescon, Float64 *qfac, Int32 *
							 kopt, Int32 *knew, Float64 *del, Float64 *step, Float64 *
							 gl, Float64 *pqw, Float64 *rstat, Float64 *w, Int32 *ifeas)
{
	/* System generated locals */
	Int32 amat_dim1, amat_offset, xpt_dim1, xpt_offset, qfac_dim1, qfac_offset, i__1, i__2;
	Float64 d__1, d__2;

	/* Builtin functions */
	double sqrt(Float64);

	/* Local variables */
	Int32 i__, j, k;
	Float64 gg, sp, ss, ww, ghg, one, sum, stp, half, vbig, bigv, vlag, ctol;
	Int32 jsav;
	Float64 temp;
	Int32 ksav;
	Float64 zero, test, vnew;
	Int32 iflag;
	Float64 vgrad, tenth, resmax, stpsav;


	/*     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the */
	/*       same as the terms with these names in SUBROUTINE LINCOB. */
	/*     KNEW is the index of the interpolation point that is going to be moved. */
	/*     DEL is the current restriction on the length of STEP, which is never */
	/*       greater than the current trust region radius DELTA. */
	/*     STEP will be set to the required step from XOPT to the new point. */
	/*     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC */
	/*       is the KNEW-th Lagrange function. It is used also for some other */
	/*       gradients of LFUNC. */
	/*     PQW provides the second derivative parameters of LFUNC. */
	/*     RSTAT and W are used for working space. Their lengths must be at least */
	/*       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the */
	/*       J-th constraint is irrelevant, active, or both inactive and relevant, */
	/*       respectively. */
	/*     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible. */

	/*     STEP is chosen to provide a relatively large value of the modulus of */
	/*       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is */
	/*       calculated too, within the trust region, that does not alter the */
	/*       residuals of the active constraints. The projected step is preferred */
	/*       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the */
	/*       original one, but the greatest violation of a linear constraint must */
	/*       be at least 0.2*DEL, in order to keep the interpolation points apart. */
	/*       The remedy when the maximum constraint violation is too small is to */
	/*       restore the original step, which is perturbed if necessary so that */
	/*       its maximum constraint violation becomes 0.2*DEL. */

	/*     Set some constants. */

	/* Parameter adjustments */
	qfac_dim1 = *n;
	qfac_offset = 1 + qfac_dim1;
	qfac -= qfac_offset;
	amat_dim1 = *n;
	amat_offset = 1 + amat_dim1;
	amat -= amat_offset;
	xpt_dim1 = *npt;
	xpt_offset = 1 + xpt_dim1;
	xpt -= xpt_offset;
	--b;
	--xopt;
	--iact;
	--rescon;
	--step;
	--gl;
	--pqw;
	--rstat;
	--w;

	/* Function Body */
	half = .5;
	one = 1.;
	tenth = .1;
	zero = 0.;
	test = *del * .2;

	/*     Replace GL by the gradient of LFUNC at the trust region centre, and */
	/*       set the elements of RSTAT. */

	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		temp = zero;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			/* L10: */
			temp += xpt[k + j * xpt_dim1] * xopt[j];
		}
		temp = pqw[k] * temp;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L20: */
			gl[i__] += temp * xpt[k + i__ * xpt_dim1];
		}
	}
	if (*m > 0) {
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
			rstat[j] = one;
			/* L30: */
			if ((d__1 = rescon[j], abs(d__1)) >= *del) {
				rstat[j] = -one;
			}
		}
		i__2 = *nact;
		for (k = 1; k <= i__2; ++k) {
			/* L40: */
			rstat[iact[k]] = zero;
		}
	}

	/*     Find the greatest modulus of LFUNC on a line through XOPT and */
	/*       another interpolation point within the trust region. */

	iflag = 0;
	vbig = zero;
	i__2 = *npt;
	for (k = 1; k <= i__2; ++k) {
		if (k == *kopt) {
			goto L60;
		}
		ss = zero;
		sp = zero;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			temp = xpt[k + i__ * xpt_dim1] - xopt[i__];
			ss += temp * temp;
			/* L50: */
			sp += gl[i__] * temp;
		}
		stp = -(*del) / sqrt(ss);
		if (k == *knew) {
			if (sp * (sp - one) < zero) {
				stp = -stp;
			}
			vlag = (d__2 = stp * sp, abs(d__2)) + stp * stp * (d__1 = sp - 
				one, abs(d__1));
		} else {
			vlag = (d__1 = stp * (one - stp) * sp, abs(d__1));
		}
		if (vlag > vbig) {
			ksav = k;
			stpsav = stp;
			vbig = vlag;
		}
L60:
		;
	}

	/*     Set STEP to the move that gives the greatest modulus calculated above. */
	/*       This move may be replaced by a steepest ascent step from XOPT. */

	gg = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
		/* Computing 2nd power */
		d__1 = gl[i__];
		gg += d__1 * d__1;
		/* L70: */
		step[i__] = stpsav * (xpt[ksav + i__ * xpt_dim1] - xopt[i__]);
	}
	vgrad = *del * sqrt(gg);
	if (vgrad <= tenth * vbig) {
		goto L220;
	}

	/*     Make the replacement if it provides a larger value of VBIG. */

	ghg = zero;
	i__2 = *npt;
	for (k = 1; k <= i__2; ++k) {
		temp = zero;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			/* L80: */
			temp += xpt[k + j * xpt_dim1] * gl[j];
		}
		/* L90: */
		ghg += pqw[k] * temp * temp;
	}
	vnew = vgrad + (d__1 = half * *del * *del * ghg / gg, abs(d__1));
	if (vnew > vbig) {
		vbig = vnew;
		stp = *del / sqrt(gg);
		if (ghg < zero) {
			stp = -stp;
		}
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L100: */
			step[i__] = stp * gl[i__];
		}
	}
	if (*nact == 0 || *nact == *n) {
		goto L220;
	}

	/*     Overwrite GL by its projection. Then set VNEW to the greatest */
	/*       value of |LFUNC| on the projected gradient from XOPT subject to */
	/*       the trust region bound. If VNEW is sufficiently large, then STEP */
	/*       may be changed to a move along the projected gradient. */

	i__2 = *n;
	for (k = *nact + 1; k <= i__2; ++k) {
		w[k] = zero;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L110: */
			w[k] += gl[i__] * qfac[i__ + k * qfac_dim1];
		}
	}
	gg = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		gl[i__] = zero;
		i__2 = *n;
		for (k = *nact + 1; k <= i__2; ++k) {
			/* L120: */
			gl[i__] += qfac[i__ + k * qfac_dim1] * w[k];
		}
		/* L130: */
		/* Computing 2nd power */
		d__1 = gl[i__];
		gg += d__1 * d__1;
	}
	vgrad = *del * sqrt(gg);
	if (vgrad <= tenth * vbig) {
		goto L220;
	}
	ghg = zero;
	i__1 = *npt;
	for (k = 1; k <= i__1; ++k) {
		temp = zero;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			/* L140: */
			temp += xpt[k + j * xpt_dim1] * gl[j];
		}
		/* L150: */
		ghg += pqw[k] * temp * temp;
	}
	vnew = vgrad + (d__1 = half * *del * *del * ghg / gg, abs(d__1));

	/*     Set W to the possible move along the projected gradient. */

	stp = *del / sqrt(gg);
	if (ghg < zero) {
		stp = -stp;
	}
	ww = zero;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = stp * gl[i__];
		/* L160: */
		/* Computing 2nd power */
		d__1 = w[i__];
		ww += d__1 * d__1;
	}

	/*     Set STEP to W if W gives a sufficiently large value of the modulus */
	/*       of the Lagrange function, and if W either preserves feasibility */
	/*       or gives a constraint violation of at least 0.2*DEL. The purpose */
	/*       of CTOL below is to provide a check on feasibility that includes */
	/*       a tolerance for contributions from computer rounding errors. */

	if (vnew / vbig >= .2) {
		*ifeas = 1;
		bigv = zero;
		j = 0;
L170:
		++j;
		if (j <= *m) {
			if (rstat[j] == one) {
				temp = -rescon[j];
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					/* L180: */
					temp += w[i__] * amat[i__ + j * amat_dim1];
				}
				bigv = max(bigv,temp);
			}
			if (bigv < test) {
				goto L170;
			}
			*ifeas = 0;
		}
		ctol = zero;
		temp = sqrt(ww) * .01;
		if (bigv > zero && bigv < temp) {
			i__1 = *nact;
			for (k = 1; k <= i__1; ++k) {
				j = iact[k];
				sum = zero;
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* L190: */
					sum += w[i__] * amat[i__ + j * amat_dim1];
				}
				/* L200: */
				/* Computing MAX */
				d__1 = ctol, d__2 = abs(sum);
				ctol = max(d__1,d__2);
			}
		}
		if (bigv <= ctol * 10. || bigv >= test) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* L210: */
				step[i__] = w[i__];
			}
			goto L260;
		}
	}

	/*     Calculate the greatest constraint violation at XOPT+STEP with STEP at */
	/*       its original value. Modify STEP if this violation is unacceptable. */

L220:
	*ifeas = 1;
	bigv = zero;
	resmax = zero;
	j = 0;
L230:
	++j;
	if (j <= *m) {
		if (rstat[j] < zero) {
			goto L230;
		}
		temp = -rescon[j];
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L240: */
			temp += step[i__] * amat[i__ + j * amat_dim1];
		}
		resmax = max(resmax,temp);
		if (temp < test) {
			if (temp <= bigv) {
				goto L230;
			}
			bigv = temp;
			jsav = j;
			*ifeas = -1;
			goto L230;
		}
		*ifeas = 0;
	}
	if (*ifeas == -1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L250: */
			step[i__] += (test - bigv) * amat[i__ + jsav * amat_dim1];
		}
		*ifeas = 0;
	}

	/*     Return the calculated STEP and the value of IFEAS. */

L260:
	return 0;
} /* qmstep_ */


};
#endif

#endif