#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(BIG-TREE-SIM, 0.0.1, mtholder@gmail.com)
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.57)


################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299 subdir-objects])

AC_CONFIG_SRCDIR([src/big_tree_sim.cpp])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX

AC_PROG_INSTALL

LIBS="$LIBS"
LDFLAGS="$LDFLAGS"
AC_SUBST(CFLAGS)
AC_SUBST(abs_top_srcdir)

AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([memmove strchr strtol])
AC_CONFIG_FILES([				\
	Makefile 					\
	src/Makefile				\
	test/Makefile
	])
AC_OUTPUT
