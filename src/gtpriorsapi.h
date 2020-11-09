SEXP pkn_ngg_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec);
SEXP pkn_ngg_rec_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec);
SEXP vnk_ngg_api_call(SEXP k, SEXP n, SEXP beta, SEXP sigma, SEXP prec);
SEXP cnk_ngg_api_call(SEXP k, SEXP n, SEXP sigma, SEXP prec);
SEXP cnk_ngg_rec_api_call(SEXP k, SEXP n, SEXP sigma, SEXP prec);

void pkn_ngg_api(int *k, int *n, double *beta, double *sigma, int *prec, double *pkn, double *radius, int *accu);
void pkn_ngg_rec_api(int *k, int *n, double *beta, double *sigma, int *prec, double *pkn, double *radius, int *accu);
void vnk_ngg_api(int *k, int *n, double *beta, double *sigma, int *prec, double *vnk, double *radius, int *accu);
void cnk_ngg_api(int *k, int *n, double *sigma, int *prec, double *cnk, double *radius, int *accu);
void cnk_ngg_rec_api(int *k, int *n, double *sigma, int *prec, double *cnk, double *radius, int *accu);

static R_NativePrimitiveArgType pkn_ngg_api_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType vnk_ngg_api_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType cnk_ngg_api_t[] = {
    INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};

static const R_CallMethodDef callMethods[] = {
	{"pkn_ngg_api_call", (DL_FUNC) &pkn_ngg_api_call, 5},
	{"pkn_ngg_rec_api_call", (DL_FUNC) &pkn_ngg_rec_api_call, 5},
	{"vnk_ngg_api_call", (DL_FUNC) &vnk_ngg_api_call, 5},
	{"cnk_ngg_api_call", (DL_FUNC) &cnk_ngg_api_call, 4},
	{"cnk_ngg_rec_api_call", (DL_FUNC) &cnk_ngg_rec_api_call, 4},
	{NULL, NULL, 0}
};

static const R_CMethodDef cMethods[] = {
	{"pkn_ngg_api", (DL_FUNC) &pkn_ngg_api, 8, pkn_ngg_api_t},
	{"pkn_ngg_rec_api", (DL_FUNC) &pkn_ngg_rec_api, 8, pkn_ngg_api_t},
	{"vnk_ngg_api", (DL_FUNC) &vnk_ngg_api, 8, vnk_ngg_api_t},
	{"cnk_ngg_api", (DL_FUNC) &cnk_ngg_api, 7, cnk_ngg_api_t},
	{"cnk_ngg_rec_api", (DL_FUNC) &cnk_ngg_rec_api, 7, cnk_ngg_api_t},
	{NULL, NULL, 0, NULL}
};

