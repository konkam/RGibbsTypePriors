#' @importFrom Rmpfr factorialMpfr igamma chooseMpfr mpfr
#' @export
#' @rdname vnk_ngg
vnk_ngg.mpfr=function(n, k, beta, sigma, prec=53){
	stopifnot(beta>=0, sigma<1, sigma>=0, all(k>=0), n>=0, prec>0, length(n)==1);

	if(class(sigma)!="mpfr"){
		sigma=mpfr(sigma, precBits=prec);
	}
	if(class(beta)!="mpfr"){
		beta=mpfr(beta, precBits=prec);
	}
	cst=exp(beta)*sigma^(k-1)/factorialMpfr(n-1, precBits=prec);
	i=0:(n-1);
	summands=chooseMpfr(n-1, i)*(-1)^i*beta^(i/sigma)*igamma(k-i/sigma, beta);
	vnk=sum(summands)*cst;
	return(list(vnk=vnk, radius=NA, accu=NA));
}

#' @importFrom Rmpfr factorialMpfr pochMpfr chooseMpfr mpfr
#' @export
#' @rdname cnk_ngg
cnk_ngg.mpfr=function(n, k, sigma, prec=53){
	stopifnot(sigma<1, sigma>=0, all(k>=0), n>=0, prec>0, length(n)==1);

	if(class(sigma)!="mpfr"){
		sigma=mpfr(sigma, precBits=prec);
	}
	cst=1/factorialMpfr(k, precBits=prec);
	i=0:k;
	summands=(-1)^i*chooseMpfr(k, i)*pochMpfr(-i*sigma, n);
	cnk=sum(summands)*cst;
	return(list(cnk=cnk, radius=NA, accu=NA));
}

#' @importFrom Rmpfr mpfr
#' @export
#' @rdname pkn_ngg
pkn_ngg.mpfr=function(k, n, beta, sigma, prec=53){
	stopifnot(beta>=0, sigma<1, sigma>=0, all(k>=0), n>=0, prec>0, length(n)==1);

	if(class(sigma)!="mpfr"){
		sigma=mpfr(sigma, precBits=prec);
	}
	if(class(beta)!="mpfr"){
		beta=mpfr(beta, precBits=prec);
	}
	pkn=sapply(k, function(i){
		vnk=vnk_ngg.mpfr(n=n, k=i, beta=beta, sigma=sigma, prec=prec)$vnk;
		cnk=cnk_ngg.mpfr(n=n, k=i, sigma=sigma, prec=prec)$cnk;
		return(vnk*cnk/sigma^i);
	});
	return(list(pkn=pkn, radius=rep(NA, times=length(k)), accu=rep(NA, times=length(k)), k=k, n=n, sigma=as.numeric(sigma), beta=as.numeric(beta)));
}


