#ifdef HAVE_ARB_H
#include <arf.h>
#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>

#define MEMO_KMAX (10000)
#define MEMO_NMAX (10000)

// Memoization structure
struct memotable_s {
	double sigma;
	double r;
	unsigned int prec;
	arb_t *table;
	int *init;
	int initialized;
};


void initialize_memoization(double sigma, double r, unsigned int prec);
void cleanup_memoization();

void vnk_ngg(arb_t *vnk, unsigned int prec, unsigned int n, unsigned int k, double beta, double sigma);
void cnk_ngg(arb_t *c, unsigned int prec, unsigned int n, unsigned int k, double sigma);
void cnk_ngg_rec(arb_t *cnk, unsigned int prec, unsigned int n, unsigned int k, double sigma);
void pkn_ngg(arb_t *p, unsigned int prec, unsigned int k, unsigned int n, double beta, double sigma);
void pkn_ngg_rec(arb_t *p, unsigned int prec, unsigned int k, unsigned int n, double beta, double sigma);
#endif
