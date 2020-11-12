#include "config.h"

#include <limits.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>

#ifdef HAVE_ARB_H
#include <mpfr.h>
#include <arb.h>
#include <arf.h>
#endif

#include "gtpriors.h"
#include "gtpriorsapi.h"


void R_init_RGibbsTypePriors(DllInfo *info){
	R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
#ifdef HAVE_ARB_H
	initialize_memoization(-1, -1, 0);
#endif
}

void R_unload_RGibbsTypePriors(DllInfo *info){
#ifdef HAVE_ARB_H
	cleanup_memoization();
#endif
}

#ifndef HAVE_ARB_H
SEXP gtprior_empty5(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec){
	return(0);
}

SEXP gtprior_empty4(SEXP k, SEXP n, SEXP sigma, SEXP prec){
	return(0);
}

#endif

#ifdef HAVE_ARB_H
void arb_to_sexp(SEXP *res, arb_t x, int prec){
	*res=PROTECT(allocVector(STRSXP, 3));
	int ndec=ceil(prec*log10(2))+50;
	char *s=malloc(ndec*sizeof(char));
	mpfr_t y;
	mpfr_init2(y, prec); 
	arf_get_mpfr(y, arb_midref(x), MPFR_RNDN);
	mpfr_sprintf(s, "%RNe", y);
	SET_STRING_ELT(*res, 0, mkChar(s));
	mpfr_clear(y);

	double radius=mag_get_d(arb_radref(x));
	sprintf(s, "%e", radius);
	SET_STRING_ELT(*res, 1, mkChar(s));

	int accu;
	if(arb_is_zero(x) || arb_rel_accuracy_bits(x)>INT_MAX){
		accu=prec;
	} else {
		accu=arb_rel_accuracy_bits(x);
	}
	sprintf(s, "%d", accu);
	SET_STRING_ELT(*res, 2, mkChar(s));

	free(s);
	UNPROTECT(1);
}

SEXP pkn_ngg_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec){
	SEXP res;
	arb_t pkn_arb;

	k=PROTECT(coerceVector(k, INTSXP));
	n=PROTECT(coerceVector(n, INTSXP));
	beta=PROTECT(coerceVector(beta, REALSXP));
	sigma=PROTECT(coerceVector(sigma, REALSXP));
	prec=PROTECT(coerceVector(prec, INTSXP));

	arb_init(pkn_arb);
	pkn_ngg(&pkn_arb, INTEGER(prec)[0], INTEGER(k)[0], INTEGER(n)[0], REAL(beta)[0], REAL(sigma)[0]);
	arb_to_sexp(&res, pkn_arb, INTEGER(prec)[0]);
	arb_clear(pkn_arb);

	UNPROTECT(5);
	return(res);
}

SEXP pkn_ngg_rec_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec){
	SEXP res;
	arb_t pkn_arb;

	k=PROTECT(coerceVector(k, INTSXP));
	n=PROTECT(coerceVector(n, INTSXP));
	beta=PROTECT(coerceVector(beta, REALSXP));
	sigma=PROTECT(coerceVector(sigma, REALSXP));
	prec=PROTECT(coerceVector(prec, INTSXP));

	arb_init(pkn_arb);
	pkn_ngg_rec(&pkn_arb, INTEGER(prec)[0], INTEGER(k)[0], INTEGER(n)[0], REAL(beta)[0], REAL(sigma)[0]);
	arb_to_sexp(&res, pkn_arb, INTEGER(prec)[0]);
	arb_clear(pkn_arb);

	UNPROTECT(5);
	return(res);
}

SEXP vnk_ngg_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec){
	SEXP res;
	arb_t vnk_arb;

	k=PROTECT(coerceVector(k, INTSXP));
	n=PROTECT(coerceVector(n, INTSXP));
	beta=PROTECT(coerceVector(beta, REALSXP));
	sigma=PROTECT(coerceVector(sigma, REALSXP));
	prec=PROTECT(coerceVector(prec, INTSXP));

	arb_init(vnk_arb);
	vnk_ngg(&vnk_arb, INTEGER(prec)[0], INTEGER(n)[0], INTEGER(k)[0], REAL(beta)[0], REAL(sigma)[0]);
	arb_to_sexp(&res, vnk_arb, INTEGER(prec)[0]);
	arb_clear(vnk_arb);

	UNPROTECT(5);
	return(res);
}

SEXP cnk_ngg_api_call(SEXP k, SEXP n, SEXP sigma, SEXP prec){
	SEXP res;
	arb_t cnk_arb;

	k=PROTECT(coerceVector(k, INTSXP));
	n=PROTECT(coerceVector(n, INTSXP));
	sigma=PROTECT(coerceVector(sigma, REALSXP));
	prec=PROTECT(coerceVector(prec, INTSXP));

	arb_init(cnk_arb);
	cnk_ngg(&cnk_arb, INTEGER(prec)[0], INTEGER(n)[0], INTEGER(k)[0],REAL(sigma)[0]);
	arb_to_sexp(&res, cnk_arb, INTEGER(prec)[0]);
	arb_clear(cnk_arb);

	UNPROTECT(5);
	return(res);
}

SEXP cnk_ngg_rec_api_call(SEXP k, SEXP n, SEXP sigma, SEXP prec){
	SEXP res;
	arb_t cnk_arb;

	k=PROTECT(coerceVector(k, INTSXP));
	n=PROTECT(coerceVector(n, INTSXP));
	sigma=PROTECT(coerceVector(sigma, REALSXP));
	prec=PROTECT(coerceVector(prec, INTSXP));

	arb_init(cnk_arb);
	cnk_ngg(&cnk_arb, INTEGER(prec)[0], INTEGER(n)[0], INTEGER(k)[0],REAL(sigma)[0]);
	arb_to_sexp(&res, cnk_arb, INTEGER(prec)[0]);
	arb_clear(cnk_arb);

	UNPROTECT(5);
	return(res);
}


void pkn_ngg_api(int *k, int *n, double *beta, double *sigma, int *prec, double *pkn, double *radius, int *accu){
	arb_t pkn_arb;

	arb_init(pkn_arb);
	pkn_ngg(&pkn_arb, prec[0], k[0], n[0], beta[0], sigma[0]);
	pkn[0]=arf_get_d(arb_midref(pkn_arb), ARF_RND_NEAR);
	radius[0]=mag_get_d(arb_radref(pkn_arb));
	if(arb_is_zero(pkn_arb)){
		accu[0]=prec[0];
	} else {
		accu[0]=arb_rel_accuracy_bits(pkn_arb);
	}

	arb_clear(pkn_arb);
}

void pkn_ngg_rec_api(int *k, int *n, double *beta, double *sigma, int *prec, double *pkn, double *radius, int *accu){
	arb_t pkn_arb;

	arb_init(pkn_arb);
	pkn_ngg_rec(&pkn_arb, prec[0], k[0], n[0], beta[0], sigma[0]);
	pkn[0]=arf_get_d(arb_midref(pkn_arb), ARF_RND_NEAR);
	radius[0]=mag_get_d(arb_radref(pkn_arb));
	if(arb_is_zero(pkn_arb)){
		accu[0]=prec[0];
	} else {
		accu[0]=arb_rel_accuracy_bits(pkn_arb);
	}

	arb_clear(pkn_arb);
}

void vnk_ngg_api(int *k, int *n, double *beta, double *sigma, int *prec, double *vnk, double *radius, int *accu){
	arb_t vnk_arb;

	arb_init(vnk_arb);
	vnk_ngg(&vnk_arb, prec[0], n[0], k[0], beta[0], sigma[0]);
	vnk[0]=arf_get_d(arb_midref(vnk_arb), ARF_RND_NEAR);
	radius[0]=mag_get_d(arb_radref(vnk_arb));
	if(arb_is_zero(vnk_arb)){
		accu[0]=prec[0];
	} else {
		accu[0]=arb_rel_accuracy_bits(vnk_arb);
	}

	arb_clear(vnk_arb);
}

void cnk_ngg_api(int *k, int *n, double *sigma, int *prec, double *cnk, double *radius, int *accu){
	arb_t cnk_arb;

	arb_init(cnk_arb);
	cnk_ngg(&cnk_arb, prec[0], n[0], k[0], sigma[0]);
	cnk[0]=arf_get_d(arb_midref(cnk_arb), ARF_RND_NEAR);
	radius[0]=mag_get_d(arb_radref(cnk_arb));
	if(arb_is_zero(cnk_arb)){
		accu[0]=prec[0];
	} else {
		accu[0]=arb_rel_accuracy_bits(cnk_arb);
	}

	arb_clear(cnk_arb);
}

void cnk_ngg_rec_api(int *k, int *n, double *sigma, int *prec, double *cnk, double *radius, int *accu){
	arb_t cnk_arb;

	arb_init(cnk_arb);
	cnk_ngg_rec(&cnk_arb, prec[0], n[0], k[0], sigma[0]);
	cnk[0]=arf_get_d(arb_midref(cnk_arb), ARF_RND_NEAR);
	radius[0]=mag_get_d(arb_radref(cnk_arb));
	if(arb_is_zero(cnk_arb)){
		accu[0]=prec[0];
	} else {
		accu[0]=arb_rel_accuracy_bits(cnk_arb);
	}

	arb_clear(cnk_arb);
}

#endif
