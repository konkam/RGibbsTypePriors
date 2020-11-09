AC_INIT

if test /usr/local != "$exec_prefix" -a /usr/local != "$prefix" ; then
  CPPFLAGS="$CFLAGS -I/usr/local/include"
  LDFLAGS="$LDFLAGS -L/usr/local/lib"
fi

AC_CONFIG_HEADER(src/config.h)

AC_CHECK_LIB(mpfr, mpfr_sprintf, [HAS_MPFR=TRUE],[])
AC_CHECK_HEADERS(mpfr.h)

AC_CHECK_LIB(flint-arb, acb_hypgeom_gamma_upper, [HAS_ARB=TRUE], [HAS_ARB=FALSE])
AC_SUBST(HAS_ARB)
AC_CHECK_HEADERS(arb.h)

PKG_CPPFLAGS="${PKG_CFLAGS} -I/usr/local/include/ "
PKG_LIBS="${PKG_LIBS} -L/usr/local/lib/ -lmpfr -lgmp -lflint-arb"
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)

AC_CONFIG_FILES([src/Makevars R/gtpriors.R])
AC_OUTPUT
