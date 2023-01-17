#include <arf.h>
#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>

void vnk_ngg(arb_t *vnk, unsigned int prec, unsigned int n, unsigned int k, double beta, double sigma){
	unsigned int i;
	int sign;
	arb_t sigma_arb, beta_arb, s, cst;
	arb_t temp0, temp1, temp2;
	arb_t bin, betapow, g, summand;
	acb_t tempC0, tempC1, tempC2;

	//init
	acb_init(tempC0);
	acb_init(tempC1);
	acb_init(tempC2);
	arb_init(sigma_arb);
	arb_init(beta_arb);
	arb_init(s);
	arb_init(temp0);
	arb_init(temp1);
	arb_init(temp2);
	arb_init(cst);
	arb_init(bin);
	arb_init(betapow);
	arb_init(g);
	arb_init(summand);

	// Sum
	arb_set_d(sigma_arb, sigma);
	arb_set_d(beta_arb, beta);

	arb_zero(*vnk);
	arb_zero(s);
	for(i=0; i<n; i++){
		sign=i%2==0?1:-1; // sign=(-1)^i

		arb_bin_uiui(bin, n-1, i, prec); // bin=C^n-1_i

		arb_set_ui(temp0, i);
		arb_div(temp1, temp0, sigma_arb, prec); // betapow=beta^(i/sigma)
		arb_pow(betapow, beta_arb, temp1, prec);

		acb_set_ui(tempC0, k);
		acb_set_ui(tempC1, i);
		acb_div_arb(tempC2, tempC1, sigma_arb, prec);
		acb_sub(tempC1, tempC0, tempC2, prec);
		acb_set_arb(tempC0, beta_arb);
		acb_hypgeom_gamma_upper(tempC2, tempC1, tempC0, 0, prec);
		acb_get_real(g, tempC2);

		arb_mul(temp0, betapow, bin, prec); // temp0=bin*betapow
		arb_mul(temp1, temp0, g, prec); // temp1=temp0*g
		arb_mul_si(summand, temp1, sign, prec); //summand=temp1*sign

		arb_set(temp0, s); // temp0=s
		arb_add(s, summand, temp0, prec); // s=temp0+summand
	}

	arb_exp(temp0, beta_arb, prec); // temp0=exp(beta)
	arb_pow_ui(temp1, sigma_arb, k-1, prec); //temp&=sigma^(k-1)
	arb_mul(temp2, temp0, temp1, prec); // temp2=temp0*temp1
	arb_fac_ui(temp1, n-1, prec); // temp1=(n-1)!
	arb_div(cst, temp2, temp1, prec); // cst=temp2/temp1;

	arb_mul(*vnk, s, cst, prec);

	// Clear
	acb_clear(tempC0);
	acb_clear(tempC1);
	acb_clear(tempC2);

	arb_clear(bin);
	arb_clear(g);
	arb_clear(betapow);
	arb_clear(cst);
	arb_clear(s);
	arb_clear(temp2);
	arb_clear(temp1);
	arb_clear(temp0);
	arb_clear(beta_arb);
	arb_clear(sigma_arb);
}

void cnk_ngg(arb_t *cnk, unsigned int prec, unsigned int n, unsigned int k, double sigma){
	int i, sign;
	arb_t s, temp0, temp1, temp2, bin, rifac, sigma_arb;
	arb_init(s);
	arb_init(bin);
	arb_init(rifac);
	arb_init(temp0);
	arb_init(temp1);
	arb_init(temp2);
	arb_init(sigma_arb);

	arb_zero(*cnk);
	arb_set_d(sigma_arb, sigma);
	for(i=0; i<k; i++){
		sign=i%2==0?1:-1;

		arb_bin_uiui(bin, k, i, prec);

		arb_mul_si(temp1, sigma_arb, -i, prec);
		arb_rising_ui(rifac, temp1, n, prec);
		
		arb_mul_si(temp0, bin, sign, prec);
		arb_mul(temp1, temp0, rifac, prec);

		arb_set(temp0, s);
		arb_add(s, temp0, temp1, prec);
	}
	arb_fac_ui(temp0, k, prec);
	arb_div(*cnk, s, temp0, prec);

	arb_clear(sigma_arb);
	arb_clear(s);
	arb_clear(bin);
	arb_clear(rifac);
	arb_clear(temp0);
	arb_clear(temp1);
	arb_clear(temp2);
}

void pkn_ngg(arb_t *p, unsigned int prec, unsigned int k, unsigned int n, double beta, double sigma){
	arb_t temp0, temp1, vnk, cnk, sigma_arb;
	if(k>n || k==0){
		arb_zero(*p);
	} else {
		arb_init(temp0);
		arb_init(temp1);
		arb_init(vnk);
		arb_init(cnk);
		arb_init(sigma_arb);

		arb_set_d(sigma_arb, sigma);
		arb_pow_ui(temp0, sigma_arb, k, prec);
		vnk_ngg(&vnk, prec, n, k, beta, sigma);
		cnk_ngg(&cnk, prec, n, k, sigma);
		arb_mul(temp1, cnk, vnk, prec);
		arb_div(*p, temp1, temp0, prec);

		arb_clear(temp0);
		arb_clear(temp1);
		arb_clear(vnk);
		arb_clear(cnk);
		arb_clear(sigma_arb);
	}
}
