#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([POY],[see VERSION],[poy4@googlegroups.com])
AC_CONFIG_SRCDIR([add.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# What are we going to substitute in the Makefile
AC_SUBST(DUSE_NATIVEDYNLINK)
AC_SUBST(DUSE_XSLT)
AC_SUBST(DUSE_LIKELIHOOD)
AC_SUBST(DUSE_PARMAP)
AC_SUBST(DUSE_CONCORDE)
AC_SUBST(CONCORDE_LIB)
AC_SUBST(CONCORDE_CFLAGS)
AC_SUBST(DUSEREADLINE)
AC_SUBST(DUSENCURSES)
AC_SUBST(DUSEPARALLEL)
AC_SUBST(DUSENOSHAREDHD)
AC_SUBST(DUSE_LARGE_MESSAGES)
AC_SUBST(DUSE_LONG_SEQUENCES)
AC_SUBST(DUSE_LARGE_ALPHABETS)
AC_SUBST(DUSE_VERIFY_COSTS)
AC_SUBST(DUSEGRAPHICS)
AC_SUBST(DUSEWIN32)
AC_SUBST(PATH_PREFIX)

AC_SUBST(OUSE_NATIVEDYNLINK)
AC_SUBST(OUSE_XSLT)
AC_SUBST(OUSE_LIKELIHOOD)
AC_SUBST(OUSEREADLINE)
AC_SUBST(OUSENCURSES)
AC_SUBST(OUSEPARALLEL)
AC_SUBST(OUSENOSHAREDHD)
AC_SUBST(OUSE_LARGE_MESSAGES)
AC_SUBST(OUSE_LONG_SEQUENCES)
AC_SUBST(OUSE_LARGE_ALPHABETS)
AC_SUBST(OUSE_VERIFY_COSTS)
AC_SUBST(OUSEGRAPHICS)
AC_SUBST(OUSEWIN32)
AC_SUBST(OUSE_SSE)
AC_SUBST(OUSE_CONCORDE)

AC_SUBST(INTERFACETOP)
AC_SUBST(MPITOP)

AC_SUBST(USEINTERFACE)
AC_SUBST(USESHAREDHD)
AC_SUBST(USEPARALLEL)
AC_SUBST(USE_LARGE_MESSAGES)
AC_SUBST(USEWIN32)
AC_SUBST(USE_LONG_SEQUENCES)
AC_SUBST(USE_LARGE_ALPHABETS)
AC_SUBST(USE_XSLT)
AC_SUBST(USE_LIKELIHOOD)
AC_SUBST(USE_PARMAP)
AC_SUBST(USE_CONCORDE)
AC_SUBST(USE_SSE)
AC_SUBST(USE_VERIFY_COSTS)
AC_SUBST(USE_UNSAFE)
AC_SUBST(USE_RTL)

AC_SUBST(POY_VERSION_NUMBER)
AC_SUBST(CC_VERSION)
AC_SUBST(OCAML)
AC_SUBST(OCAMLC)
AC_SUBST(OCAMLOPT)
AC_SUBST(OCAMLDEP)
AC_SUBST(CAMLP4ORF)
AC_SUBST(OCAMLTOP)
AC_SUBST(OCAMLLEX)
AC_SUBST(OCAMLYACC)
AC_SUBST(OCAMLBUILD)
AC_SUBST(MERCURIAL)
AC_SUBST(LIBS)
AC_SUBST(OCAMLFLAGS)
AC_SUBST(EXTRAS)
AC_SUBST(NATIVEDYNLINK)
AC_SUBST(USE_NATIVEDYNLINK)
AC_SUBST(OCAMLVERSION)
AC_SUBST(CONFIG_DATE)

AC_SUBST(CC_NOASNEEDED)
AC_SUBST(PROFILE_OCAML)

CONFIG_DATE=`date`
USEWIN32="false"
OS_INFORMATION=`../gnu/config.guess`

# Set the default; OSX adds usr later to avoid conflict with sysroot/isysroot
case ${OS_INFORMATION} in
    *-apple-*)
        PATH_PREFIX="/"
        USEINTERFACE="ncurses"
        ;;
    *-cygwin*)
        USEINTERFACE="pdcurses"
        PATH_PREFIX="/usr"
        ;;
    *)
        USEINTERFACE="ncurses"
        PATH_PREFIX="/usr"
        ;;
esac

AC_ARG_ENABLE([interface], [AS_HELP_STRING([--enable-interface],
              [Force POYs user interface. Valid options are flat, ncurses, pdcurses, readline, and html. The default is ncurses, if ncurses is not found, then the readline interface is used, if readline is not present either, then choose the flat interface.])],
              [USEINTERFACE=$enableval;INTERFACE_SELECTED="true"],
              [INTERFACE_SELECTED="false"; DUSENCURSES="-DUSENCURSES"; OUSENCURSES=",use_ncurses"])

AC_ARG_WITH([extras], [AS_HELP_STRING([--with-extras], 
            [Add the following object files to the final linking rules.])],
            [EXTRAS=$withval],[])

AC_ARG_WITH(version-number,
            [AS_HELP_STRING([--with-version-number], [Use the configuration provided version number (used for release purposes only).])],
            [MERCURIAL="echo $withval"],
            [AC_CHECK_PROG([MERCURIAL],[hg],[`hg id -i`],[`cat ./VERSION | awk '{print $NF}'`])])

AC_ARG_WITH([ocaml-flags], 
            [AS_HELP_STRING([--with-ocaml-flags], [Add the following flags to the ocaml compiler.])],
            [OCAMLFLAGS=$withval],[])

AC_ARG_ENABLE([rtlout], [AS_HELP_STRING([--enable-rtlout],
              [Used to generate RTL output that can be used for egypt to create graphviz of C routines (in _build directory). (debug use only).])],
              [if test x"$enableval" = x"no"; then
                    USE_RTL="false"
                else
                    USE_RTL="true"
                fi],
              [USE_RTL="false"])

AC_ARG_ENABLE([unsafe], [AS_HELP_STRING([--enable-unsafe],
              [When set, we increase inline level, compile with unsafe, and turn off assertions. Used for final releases to increase speed.])],
              [if test x"$enableval" = x"no"; then
                    USE_UNSAFE="false"
                else
                    USE_UNSAFE="true"
                fi],
              [USE_UNSAFE="false"])

AC_ARG_ENABLE([mpi], [AS_HELP_STRING([--enable-mpi], 
              [Turns on parallel support using MPI. When using MPI remember to set CC to the appropriate mpicc compiler.])], 
              [USEPARALLEL="true"; DUSEPARALLEL="-DUSEPARALLEL"; OUSEPARALLEL=",use_parallel"; MPITOP="Mpi"], 
              [USEPARALLEL="false"])

AC_ARG_ENABLE([large-messages], [AS_HELP_STRING([--enable-large-messages], 
              [When --enable-mpi is selected, this option turns on large messages for tree exchanges. Large messages reduce considerably computation of trees when different processes exchange them, but increase tremendously the message size. Use this option when message size (and the memory consumed to send and receive it) is not a possible bottleneck.])], 
              [USE_LARGE_MESSAGES="true"; DUSE_LARGE_MESSAGES="-DUSE_LARGE_MESSAGES"; OUSE_LARGE_MESSAGES=",use_large_messages"], 
              [USE_LARGE_MESSAGES="false"])

AC_ARG_ENABLE([parmap], [AS_HELP_STRING([--enable-parmap], 
              [Enables parmap through ocamlfind tree level parallelism in building and swappingr.])], 
              [USE_PARMAP="true";DUSE_PARMAP="-DUSE_PARMAP"], 
              [USE_PARMAP="false"])

AC_ARG_ENABLE([shared-hd], [AS_HELP_STRING([--enable-shared-hd],
              [Force POY to support live broadcast of file contents when running in parallel. The default is false.])],
              [USESHAREDHD="true"; DUSENOSHAREDHD="-DUSENOSHAREDHD"; OUSENOSHAREDHD=",use_nosharedhd"],
              [USESHAREDHD="false"])

AC_ARG_ENABLE([long-sequences], [AS_HELP_STRING([--enable-long-sequences], 
              [Enable the analysis of sequences longer than 16383 elements (The limit would be  2^31 or 2^63 for 32 bit and 64 bit architectures respectively). When using this option, more memory will be needed for the default search strategies. If you are analyzing fragments shorted than, say 12.000 bp, you are unlikely to need this option.])], 
              [USE_LONG_SEQUENCES="true"; DUSE_LONG_SEQUENCES="-DUSE_LONG_SEQUENCES"; OUSE_LONG_SEQUENCES=",use_long_sequences"], 
              [USE_LONG_SEQUENCES="false"])

AC_ARG_ENABLE([profiler], [AS_HELP_STRING([--enable-profiler],
              [Due to a bug in ocamlbuild, we must pass the -p option ourselves to allow profiling with gprof])],
              [if test x"$enableval" = x"no"; then
                    PROFILE_OCAML="false"
                else
                    PROFILE_OCAML="true"
                fi],
              [PROFILE_OCAML="false"])

AC_ARG_ENABLE([xslt], [AS_HELP_STRING([--enable-xslt],
              [Enable the postprocessing of XML output using XSLT templates.  This option requires libxml2 and libxslt.])],
              [if test x"$enableval" = x"no"; then
                    USE_XSLT="false"
                else
                    USE_XSLT="true"
                fi],
              [USE_XSLT="false"])

AC_ARG_ENABLE([likelihood], [AS_HELP_STRING([--enable-likelihood],
              [Enable/Disable the use of static and dynamic likelihood functions. We will automatically detect and enable likelihood if we can --we require lapack routines to be installed.])],
              [if test x"$enableval" = x"no"; then
                    FORCE_LIKELIHOOD="no"; USE_LIKELIHOOD="false"
               else
                    FORCE_LIKELIHOOD="yes"; USE_LIKELIHOOD="true"
               fi],
               [USE_LIKELIHOOD="true"; FORCE_LIKELIHOOD="no"])

#AC_ARG_ENABLE([concorde], [AS_HELP_STRING([--enable-concorde],
#              [Enable the use of "CHLINKERN" and "GREEDYLK" tsp median solver. This
#               option requires concorde package.])],
#                [USE_CONCORDE="true"; DUSE_CONCORDE="-DUSE_CONCORDE";
#                 OUSE_CONCORDE=",use_concorde"], [USE_CONCORDE="false"])

AC_ARG_ENABLE([sse], [AS_HELP_STRING([--enable-sse],
              [Enable/Disable the Streaming SIMD Extensions for likelihood characters.
               On by default and auto-selects SSE3 or SSE4.])],
               [ if test x"$enableval" = x"no"; then
                    USE_SSE="false"
                 else
                    USE_SSE="true"
                 fi],
                 [USE_SSE="false"])

AC_ARG_ENABLE([large-alphabets], [AS_HELP_STRING([--enable-large-alphabets],
              [Enable the analysis of alphabets larger than 255 elements. When using this option, more memory will be needed to store each tree (4 or 8 times more for 32 and 64 bit computers respectively). This option is only needed for large alphabets if you use custom alphabets. Normal nucleotide, aminoacid, or morphological analyses are nor affected by this option.])],
              [USE_LARGE_ALPHABETS="true"; DUSE_LARGE_ALPHABETS="-DUSE_LARGE_ALPHABETS"; OUSE_LARGE_ALPHABETS=",use_large_alphabets"],
              [USE_LARGE_ALPHABETS="false"])

AC_ARG_ENABLE([cost-verification], [AS_HELP_STRING([--enable-cost-verification],
              [Check every pairwise sequence alignment cost. This is only intended for development and debugging purposes.])],
              [USE_VERIFY_COSTS="true"; DUSE_VERIFY_COSTS="-DUSE_VERIFY_COSTS"; OUSE_VERIFY_COSTS=",use_verify_costs"],
              [USE_VERIFY_COSTS="false"])

AC_ARG_WITH([concorde_dir], [AS_HELP_STRING([--with-concorde-dir=<DIR>],
            [set the path to concorde library])],
            [ if test x"$withval" != x"no"; then
                concorde_dir=$withval
              fi],[])

AC_ARG_WITH([prefix], [AS_HELP_STRING([--with-prefix=<DIR>],
            [Set the build path for libraries; it should contain <DIR>/include and <DIR>/lib (for windows/linux), and in OSX is equivlent to setting isysroot (requiring <DIR>/usr to exist).])],
            [ PATH_PREFIX=$withval])

AC_ARG_WITH([xslt-config], [AS_HELP_STRING([--with-xslt-config=<PROG>],
            [Select a different xslt-config program to detect the configuration options of libxslt.])],
            [XSLT_CONFIG=$withval],[XSLT_CONFIG="xslt-config"])



# System dependent options -----
case ${OS_INFORMATION} in
    x86_64*)
        LIBS="-lm $LIBS"
        OCAMLFLAGS="-fno-PIC $OCAMLFLAGS" 
        CFLAGS="-fPIC $CFLAGS"
        ;;
    *-apple-*)
        LIBS="-isysroot $PATH_PREFIX -lm $LIBS"
        CFLAGS="-I $PATH_PREFIX/usr/include/malloc -isysroot $PATH_PREFIX $CFLAGS"
        PATH_PREFIX="$PATH_PREFIX/usr"
        ;;
    *-cygwin*)
        # m(ath) is located in msvcrt
        USEWIN32="true"
        DUSEWIN32="-DUSEWIN32"
        OUSEWIN32=",use_win32"
        CFLAGS="-mwin32 $CFLAGS"
        # so AC_FUNC_MALLOC does not redine malloc to rpl_malloc
        ac_cv_func_realloc_0_nonnull=yes
        ac_cv_func_malloc_0_nonnull=yes
        # WIN doesn't use mpicc to build/include libs; assume they're in the path
        AS_IF([test $USEPARALLEL = "true"],[LIBS="-lmpi $LIBS"])
        ;;
    *)
        LIBS="-lm $LIBS"
        ;;
esac


# This is a bugfix for OCAMLBUILD; in recent versions of gcc the order of the
# link arguments matters --the default changes from below to include libraries
# as needed. OCAMLBUILD just throws them in the constructed command.
#       Bug reported : http://caml.inria.fr/mantis/view.php?id=5441
case ${OS_INFORMATION} in
    *-cygwin*)
        AC_MSG_CHECKING([whether we need to add --no-as-needed linking option])
        CC_NOASNEEDED="false"
        AC_MSG_RESULT([not needed for flexlink]);
        ;;
    *)
        oCFLAGS="$CFLAGS"
        CFLAGS="$FLAGS -Wl,--no-as-needed"
        AC_MSG_CHECKING([whether we need to add --no-as-needed linking option])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],[[]])], 
                          [AC_MSG_RESULT([yes]); CC_NOASNEEDED="true"],
                          [AC_MSG_RESULT([no]); CC_NOASNEEDED="false"]) 
        CFLAGS=$oCFLAGS
        ;;
esac

AC_MSG_CHECKING([Whether to enable Concorde])
if test -f "$concorde_dir/concorde.a"; then
   AC_MSG_RESULT([yes])
   CONCORDE_LIB="$concorde_dir/concorde.a"
   CONCORDE_CFLAGS="-I$concorde_dir/INCLUDE"
   CFLAGS="$CFLAGS $CONCORDE_CFLAGS -DUSE_CONCORDE"
   LIBS="-L$concorde_dir $LIBS"
   USE_CONCORDE="true"
   OUSE_CONCORDE=",use_concorde"
   DUSE_CONCORDE="-DUSE_CONCORDE"
else
   AC_MSG_RESULT([no])
   USE_CONCORDE="false"
   CONCORDE_LIB=""
   CONCORDE_CFLAGS=""
   OUSE_CONCORDE=""
fi

USE_SSE="false"     #DISABLE ALL SSE until it's faster than normal.
USED_SSE="false"
AS_IF([test $USE_SSE = "true"],
    [AC_MSG_CHECKING([SSE4 support]); CFLAGS_TMP="$CFLAGS"; CFLAGS="-msse4";
     AC_RUN_IFELSE([AC_LANG_PROGRAM(
                [[
                    #include <smmintrin.h>
                ]],
                [[ 
                    __m128d res = _mm_dp_pd( _mm_set_pd(1.0,2.0), _mm_set_pd(2.0,1.0), 0x31);
                ]])],
            [CFLAGS="-msse4 $CFLAGS_TMP";USED_SSE="true";OUSE_SSE=",use_sse4";AC_MSG_RESULT([yes])],
            [CFLAGS="$CFLAGS_TMP";AC_MSG_RESULT([no])],
            [CFLAGS="$CFLAGS_TMP";AC_MSG_RESULT([disabled, cross-compiling])])],
    [AC_MSG_CHECKING([SSE4 support]);AC_MSG_RESULT([disabled])])

if test x"$USED_SSE" = x"false"; then
    AS_IF([test $USE_SSE = "true"],
        [AC_MSG_CHECKING([SSE3 support]); CFLAGS_TMP="$CFLAGS"; CFLAGS="-msse3";
         AC_RUN_IFELSE([AC_LANG_PROGRAM(
                    [[
                        #include <pmmintrin.h>
                    ]], 
                    [[ 
                        __m128d res = _mm_hadd_pd( _mm_set_pd(1.0,2.0), _mm_set_pd(1.0,2.0));
                    ]])],
                [CFLAGS="-msse3 $CFLAGS_TMP";OUSE_SSE=",use_sse3";AC_MSG_RESULT([yes])],
                [CFLAGS="$CFLAGS_TMP";AC_MSG_RESULT([no])],
                [CFLAGS="$CFLAGS_TMP";AC_MSG_RESULT([disabled, cross-compiling])])],
        [AC_MSG_CHECKING([SSE3 support]);AC_MSG_RESULT([disabled])])
fi

# Time to verify if we can compile with xslt
# Use XSLT-Config if it exists; Otherwise, assume it's in the path
AC_CHECK_PROG([HAVE_XSLTCONFIG],[$XSLT_CONFIG],yes,no)
if test "$USE_XSLT" = "true"; then
    oLIBS="$LIBS"
    if test "$HAVE_XSLTCONFIG" = "yes"; then
        LIBS="`$XSLT_CONFIG --libs` $LIBS"
        CFLAGS="`$XSLT_CONFIG --cflags` $CFLAGS"
    else
        LIBS="-lxslt -lxml2 $LIBS"
    fi

    AC_CHECK_LIB([xml2],[main],
        [AC_CHECK_LIB([xslt],[main],
            [AC_DEFINE(HAVE_XSLT)
             CFLAGS="$CFLAGS $XSLTFLAGS"
             DUSE_XSLT="-DUSE_XSLT";OUSE_XSLT=",use_xslt"],
            [LIBS="$oLIBS"])],
        [LIBS="$oLIBS"])
fi

AS_IF( [test $USEPARALLEL = "true"],
       AS_IF( [test $USEINTERFACE = "flat"], [],
              AC_MSG_RESULT([Setting interface to flat for mpi compilation]);USEINTERFACE="flat"))

AS_IF( [test $USE_PARMAP = "true"],
       AS_IF( [test $USEINTERFACE = "ncurses"],
              AC_MSG_RESULT([Setting interface to readline for parmap compilation]);USEINTERFACE="readline"))

# Checking for libraries according to the user's options.
AS_IF([test $USEINTERFACE = "ncurses"],
      [AC_CHECK_LIB([ncurses], [main],
            [LIBS="-lncurses $LIBS";DUSENCURSES="-DUSENCURSES";INTERFACETOP="NcursesML";OUSENCURSES=",use_ncurses";AC_DEFINE(HAVE_LIBNCURSES)],
            [AS_IF([test $INTERFACE_SELECTED = "true"],
                [AC_MSG_FAILURE([The ncurses interface requires the ncurses library. See ./configure --help for more information.])],
                [USEINTERFACE="readline"])])],
      [])

AS_IF([test $USEINTERFACE = "pdcurses"],
      [AC_CHECK_LIB([pdcurses], [main],
            [LIBS="-lpdcurses $LIBS";USEINTERFACE="ncurses";DUSENCURSES="-DUSENCURSES";INTERFACETOP="NcursesML";OUSENCURSES=",use_ncurses";AC_DEFINE(HAVE_LIBNCURSES)],
            [AS_IF([test $INTERFACE_SELECTED = "true"],
                [AC_MSG_FAILURE([The pdcurses interface requires the pdcurses library. See ./configure --help for more information.])],
                [USEINTERFACE="readline"])])],
      [])

AS_IF([test $USEINTERFACE = "readline"],
        [AC_CHECK_LIB([readline], [readline],
            [LIBS="-lreadline $LIBS";AC_DEFINE(HAVE_LIBREADLINE)
             USEINTERFACE="readline";DUSEREADLINE="-DUSEREADLINE";OUSEREADLINE=",use_readline"],
            [AS_IF([test $INTERFACE_SELECTED = "true"],
                [AC_MSG_FAILURE([The readline interface requires the readline library. See ./configure --help for more information.])],
                [USEINTERFACE="flat"])])],
        [])

# Check for malloc, if not there, this is most likely a Mac OS X.
AC_CHECK_HEADER(malloc.h,[],
                [AC_CHECK_HEADER("$PATH_PREFIX/include/malloc/malloc.h",
                                [CFLAGS="-I $PATH_PREFIX/include/malloc/ $CFLAGS"],
                                [AC_MSG_FAILURE([Could not find malloc.h])])])

AC_CHECK_PROGS([OCAML],[ocaml])
AC_CHECK_PROGS([OCAMLC], [ocamlc.opt ocamlc])
AC_CHECK_PROGS([OCAMLYACC], [ocamlyacc.opt ocamlyacc])
AC_CHECK_PROGS([OCAMLLEX], [ocamllex.opt ocamllex])
AC_CHECK_PROGS([OCAMLOPT], [ocamlopt.opt ocamlopt])
AC_CHECK_PROGS([OCAMLDEP], [ocamldep.opt ocamldep])
AC_CHECK_PROGS([OCAMLBUILD], [ocamlbuild.native ocamlbuild])
AC_CHECK_PROGS([CAMLP4ORF], [camlp4orf])
AC_CHECK_PROGS([OCAMLTOP], [ocamlmktop])

# Check the version of OCaml
AC_MSG_CHECKING([for OCaml version])
OCAMLVERSION=`$OCAMLOPT -version`
AC_MSG_RESULT([$OCAMLVERSION])
case $OCAMLVERSION in 
    3.10.2*)
        case ${OS_INFORMATION} in
            *-cygwin*)
                AC_MSG_RESULT([no])
                AC_MSG_FAILURE([Compiling POY in windows requires OCaml version 3.11.0 or later. Please update OCaml and configure again.])
                ;;
            *)
                AC_MSG_RESULT([yes])
                ;;
        esac
        ;;
    3.\(0*\)|\(10\).*)
        AC_MSG_RESULT([no])
        AC_MSG_FAILURE([POY requires OCaml version 3.10.2 or later to compile.  Please upgrade to the latest version. (You can get it from http://www.ocaml.org)])
        ;;
    *)
        NATIVEDYNLINK="dynlink.cmxa"
        AC_MSG_CHECKING([if OCaml supports dynamic linking in this architecture])
        if ${OCAMLOPT} -shared -o str.cmxs str.cmxa ; then
            AC_MSG_RESULT([yes])
            USE_NATIVEDYNLINK="true"
            OUSE_NATIVEDYNLINK=",use_nativedynlink"
            DUSE_NATIVEDYNLINK="-DUSE_NATIVEDYNLINK"
        else 
            AC_MSG_RESULT([no])
            USE_NATIVEDYNLINK="false"
        fi
	    rm -f str.cmxs*
        ;;
esac

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h float.h limits.h stddef.h stdint.h stdlib.h string.h strings.h sys/ioctl.h sys/time.h unistd.h],[],[],[AC_INCLUDES_DEFAULT])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

AC_PROG_GCC_TRADITIONAL
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT

AC_SEARCH_LIBS([zlibVersion],[z],,[AC_MSG_FAILURE( [I cannot find libz] )])

#this will only pass if m is already included; we include from OS switch above,
AC_CHECK_FUNCS([floor gettimeofday memset pow sqrt sin exp])

# Time to handle the lapack checks; first we check if we need to check, then we
# check the fortran routines, then lapack, then fail if there is an
# inconsistency (the case being if we explicitly enable and it fails).
AC_MSG_CHECKING([to enable likelihood])
case ${OS_INFORMATION} in
    *-apple-*)
        AS_IF([test $USE_LIKELIHOOD = "true"],
                [LIBS="-framework veclib $LIBS";USE_LIKELIHOOD="true";DUSE_LIKELIHOOD="-DUSE_LIKELIHOOD";
                 OUSE_LIKELIHOOD=",use_likelihood";AC_MSG_RESULT([enabled]);],
                [AC_MSG_RESULT([disabled])])
        ;;
    *)
        AS_IF([test $USE_LIKELIHOOD = "false"],
            [AC_MSG_RESULT([disabled])],
            [AC_MSG_RESULT([attempting])
             AC_CHECK_LIB(gfortran,main,LFORTRAN="-lgfortran",
                AC_CHECK_LIB(g2c,main,LFORTRAN="-lg2c",
                  AC_CHECK_LIB(f2c,main,LFORTRAN="-lf2c",
                    [USE_LIKELIHOOD="false"])))])

        AS_IF([test $USE_LIKELIHOOD = "false"],
            [],
            [AC_CHECK_LIB([lapack],[main],
                [USE_LIKELIHOOD="true";DUSE_LIKELIHOOD="-DUSE_LIKELIHOOD";OUSE_LIKELIHOOD=",use_likelihood";
                 LIBS="-llapack -lblas $LFORTRAN $LIBS"], [])])

        #if we cannot enable likelihood; check if this was forced!
        AC_MSG_CHECKING([if we enabled likelihood analysis])
        if [(test $FORCE_LIKELIHOOD = "yes") && (test $USE_LIKELIHOOD = "false")]; then
            AC_MSG_RESULT( [no] )
            AC_MSG_FAILURE( [I cannot find LAPACK or FORTRAN routines to enable a likelihood analysis.] )
        elif [test $USE_LIKELIHOOD = "true"]; then
            AC_MSG_RESULT( [yes] )
        else
            AC_MSG_RESULT( [disabled] )
        fi
        ;;
esac

CC_VERSION="`$CC --version`"
POY_VERSION_NUMBER="$MERCURIAL"

AC_CONFIG_FILES([Makefile _tags myocamlbuild.ml buildNumber.ml compileFlags.ml poy.mltop])
AC_OUTPUT
