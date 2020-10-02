#include <stdio.h>
#include <arb.h>
#include <arf.h>
#include <R_ext/Utils.h>

#include "gtpriors.h"

//void R_CheckUserInterrupt(void);

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
