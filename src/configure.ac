#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([POY 4],[see VERSION],[poy4@googlegroups.com])
AC_CONFIG_SRCDIR([add.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# What are we going to substitute in the Makefile
AC_SUBST(DUSE_NATIVEDYNLINK)
AC_SUBST(DUSE_XSLT)
AC_SUBST(DUSE_LIKELIHOOD)
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
AC_SUBST(OUSE_SSE3)
AC_SUBST(OUSE_CONCORDE)

AC_SUBST(USEINTERFACE)
AC_SUBST(USESHAREDHD)
AC_SUBST(USEPARALLEL)
AC_SUBST(USE_LARGE_MESSAGES)
AC_SUBST(USEWIN32)
AC_SUBST(USE_LONG_SEQUENCES)
AC_SUBST(USE_LARGE_ALPHABETS)
AC_SUBST(USE_XSLT)
AC_SUBST(USE_LIKELIHOOD)
AC_SUBST(USE_CONCORDE)
AC_SUBST(USE_SSE3)
AC_SUBST(USE_VERIFY_COSTS)

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

CONFIG_DATE=`date`
USEWIN32="false"
OS_INFORMATION=`../gnu/config.guess`
USE_LIKELIHOOD="true"

# System dependent options
case ${OS_INFORMATION} in
    x86_64*)
        LIBS="-lm $LIBS"
        OCAMLFLAGS="-fno-PIC $OCAMLFLAGS" 
        CFLAGS="-fPIC $CFLAGS"
        ;;
    *-apple-*)
        LIBS="-lm $LIBS"
        CFLAGS="-I /usr/include/malloc $CFLAGS"
        ;;
    *-cygwin*)
        USEWIN32="true"
        DUSEWIN32="-DUSEWIN32"
        OUSEWIN32=",use_win32"
        CFLAGS="-mwin32 -mno-cygwin $CFLAGS"
        # m(ath) is located in msvcrt
        LIBS="$LIBS"
        # so AC_FUNC_MALLOC does not redine malloc to rpl_malloc
        ac_cv_func_malloc_0_nonnull=no
        ac_cv_func_realloc_0_nonnull=no
        ;;
    *)
        ;;
esac

INTERFACE_SELECTED="true"

AC_ARG_ENABLE([interface], [AS_HELP_STRING([--enable-interface],
              [Force POY's user interface. Valid options are flat, ncurses, readline, and html. The default is ncurses, if ncurses is not found, then the readline interface is used, if readline is not present either, then choose the flat interface.])],
              [USEINTERFACE=$enableval],
              [USEINTERFACE="ncurses"; INTERFACE_SELECTED="false"; DUSENCURSES="-DUSENCURSES"; OUSENCURSES=",use_ncurses"])

AC_ARG_WITH(extras, [AS_HELP_STRING([--with-extras], 
            [Add the following object files to the final linking rules (needed for windows).])],
            [EXTRAS=$withval],[])

AC_CHECK_PROGS([MERCURIAL], [hg])

AC_ARG_WITH(version-number, [AS_HELP_STRING([--with-version-number], [Use the configuration provided version number (used for release purposes only).])], [MERCURIAL="echo $withval"],[MERCURIAL="$MERCURIAL id -i"])

AC_ARG_WITH(ocaml-flags, [AS_HELP_STRING([--with-ocaml-flags], 
            [Add the following flags to the ocaml compiler.])],
            [OCAMLFLAGS=$withval],[])

AC_ARG_ENABLE([mpi], [AS_HELP_STRING([--enable-mpi], 
              [Turns on parallel support using MPI. When using MPI remember to set CC to the appropriate mpicc compiler.])], 
              [USEPARALLEL="true"; DUSEPARALLEL="-DUSEPARALLEL"; OUSEPARALLEL=",use_parallel"], 
              [USEPARALLEL="false"])

AC_ARG_ENABLE([large-messages], [AS_HELP_STRING([--enable-large-messages], 
              [When --enable-mpi is selected, this option turns on large messages for tree exchanges. Large messages reduce considerably computation of trees when different processes exchange them, but increase tremendously the message size. Use this option when message size (and the memory consumed to send and receive it) is not a possible bottleneck.])], 
              [USE_LARGE_MESSAGES="true"; DUSE_LARGE_MESSAGES="-DUSE_LARGE_MESSAGES"; OUSE_LARGE_MESSAGES=",use_large_messages"], 
              [USE_LARGE_MESSAGES="false"])

AC_ARG_ENABLE([shared-hd], [AS_HELP_STRING([--enable-shared-hd],
              [Force POY to support live broadcast of file contents when running in parallel. The default is false.])], [USESHAREDHD="true"; DUSENOSHAREDHD="-DUSENOSHAREDHD"; OUSENOSHAREDHD=",use_nosharedhd"], [USESHAREDHD="false"])

AC_ARG_ENABLE([long-sequences], [AS_HELP_STRING([--enable-long-sequences], 
              [Enable the analysis of sequences longer than 16383 elements (The limit would be  2^31 or 2^63 for 32 bit and 64 bit architectures respectively). When using this option, more memory will be needed for the default search strategies. If you are analyzing fragments shorted than, say 12.000 bp, you are unlikely to need this option.])], 
              [USE_LONG_SEQUENCES="true"; DUSE_LONG_SEQUENCES="-DUSE_LONG_SEQUENCES"; OUSE_LONG_SEQUENCES=",use_long_sequences"], 
              [USE_LONG_SEQUENCES="false"])

AC_ARG_ENABLE([xslt], [AS_HELP_STRING([--enable-xslt],
              [Enable the postprocessing of XML output using XSLT templates.  This option requires libxml2.])],
               [USE_XSLT="true"; DUSE_XSLT="-DUSE_XSLT"; OUSE_XSLT=",use_xslt"],
               [USE_XSLT="false"])

AC_ARG_ENABLE([likelihood], [AS_HELP_STRING([--disable-likelihood],
              [Disable the use of static and dynamic likelihood functions. We
               will automatically detect and enable likelihood if we can --we
               require lapack routines to be installed.])],
                [USE_LIKELIHOOD="false"])

#AC_ARG_ENABLE([concorde], [AS_HELP_STRING([--enable-concorde],
#              [Enable the use of "CHLINKERN" and "GREEDYLK" tsp median solver. This
#               option requires concorde package.])],
#                [USE_CONCORDE="true"; DUSE_CONCORDE="-DUSE_CONCORDE";
#                 OUSE_CONCORDE=",use_concorde"], [USE_CONCORDE="false"])

#AC_ARG_ENABLE([sse3], [AS_HELP_STRING([--enable-sse3],
#              [Enable the Streaming SIMD Extensions for likelihood characters.])],
#               [USE_SSE3="true"],[USE_SSE3="false"])

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


AC_MSG_CHECKING([for concorde.a])
if test -f "$concorde_dir/concorde.a"; then
   AC_MSG_RESULT([yes])
   AC_MSG_RESULT([Using Concorde])
   CONCORDE_LIB="$concorde_dir/concorde.a"
   CONCORDE_CFLAGS="-I$concorde_dir/INCLUDE"
   CFLAGS="$CFLAGS $CONCORDE_CFLAGS -DUSE_CONCORDE"
   LIBS="-L$concorde_dir $LIBS"
   OUSE_CONCORDE=",use_concorde"
   DUSE_CONCORDE="-DUSE_CONCORDE"
else
   AC_MSG_RESULT([no])
   CONCORDE_LIB=""
   CONCORDE_CFLAGS=""
   OUSE_CONCORDE=""
   AC_MSG_RESULT( [Not Using Concorde ] )
fi


#AS_IF([test $USE_CONCORDE = "true"],
#				[
#                CFLAGS="$CFLAGS $CONCORDE_CFLAGS -DUSE_CONCORDE";
#                CFLAGS="$CFLAGS -I/$HOME/concorde/INCLUDE/ -DUSE_CONCORDE";
#                LIBS="$HOME/concorde/concorde.a $LIBS";
#                LIBS="$CONCORDE_LIB $LIBS";
#                OUSE_CONCORDE=",use_concorde";
#                DUSE_CONCORDE="-DUSE_CONCORDE";
#               ],
#                 [AC_MSG_RESULT( [Not Using Concorde ])] )

# check all the individual functions for the lapack routines
#       AS_IF([test $USE_LIKELIHOOD = "true"],
#           [AC_CHECK_LIB([lapack],[dsyev], [AC_CHECK_LIB([lapack],[dgeev], [AC_CHECK_LIB([lapack],[dgemm],
#           [AC_CHECK_LIB([lapack],[dgetri], [AC_CHECK_LIB([lapack],[dgetrf],
#                   [AC_DEFINE(HAVE_LIBLAPACK,1,"lapack confirmation") LIBS=[-llapack $LIBS] ],
#                   [AC_MSG_FAILURE( [Missing lapack routine: dgetrf_.])])], [AC_MSG_FAILURE( [Missing lapack routine: dgetri_.])])],
#                   [AC_MSG_FAILURE( [Missing lapack routine: dgemm_.])])], [AC_MSG_FAILURE( [Missing lapack routine: dgeev_.])])],
#                   [AC_MSG_FAILURE( [Missing lapack routine: dsyev_.])])],
#           [AC_MSG_RESULT([Not including LAPACK])])

#case ${OS_INFORMATION} in
#    *-apple-*)
#        AS_IF([test $USE_SSE3 = "true"],
#            [AC_CHECK_HEADER([pmmintrin.h],
#                [CFLAGS="-msse3 $CFLAGS";OUSE_SSE3=",use_sse3"],
#                [AC_MSG_FAILURE([Cannot compile with SSE3 instructions])])],
#            AC_MSG_RESULT([not using SSE3]))
#        ;;
#    *)
#        AS_IF([test $USE_SSE3 = "true"],
#            [AC_CHECK_HEADER([pmmintrin.h],
#                   [AC_CHECK_FUNC([posix_memalign],
#                       [CFLAGS="-msse3 $CFLAGS";OUSE_SSE3=",use_sse3"],
#                     [AC_MSG_FAILURE([posix_memalign not found])])],
#                 [AC_MSG_FAILURE([Cannot compile with SSE3 instructions])])],
#            AC_MSG_RESULT([Not using SSE3]))
#        ;;
#esac

# Check if the user wants SSE3 instructions, copy cflags to temporary variable
# to compile test SSE3
#CFLAGS_TMP=$CFLAGS
#CFLAGS="-msse3"
#AS_IF([test $USE_SSE3 = "true"],
#    [AC_MSG_CHECKING([SSE3 support...]);
#     AC_COMPILE_IFELSE([
#        #include <pmmintrin.h>
#        __m128d _add(__m128d v) { return _mm_hadd_pd( v, v ); }], 
#        [AC_CHECK_FUNC(posix_memalign, [CFLAGS="-msse3 $CFLAGS_TMP";OUSE_SSE3=",use_sse3";AC_MSG_RESULT([yes]) ], 
#                                       [AC_MSG_FAILURE([posix_memalign memory alloctaion not found])])],
#        [AC_MSG_FAILURE([Cannot compile SSE3 instructions])]
#        )]
#    )

# Checking for libraries according to the user's options.
AS_IF([test $USEINTERFACE = "ncurses"], 
      [AC_CHECK_LIB([ncurses], [scrollok], [LIBS="-lncurses $LIBS";DUSENCURSES="-DUSENCURSES";OUSENCURSES=",use_ncurses";AC_DEFINE(HAVE_LIBNCURSES)], 
      [AS_IF([test $INTERFACE_SELECTED = "true"], 
             [AC_MSG_FAILURE([The ncurses interface requires the ncurses library. Either install it, or configure with --enable-interface=readline interface. See ./configure --help for more information.])], 
             [USEINTERFACE="readline"; DUSEREADLINE="-DUSEREADLINE"; OUSEREADLINE=",use_readline"])])],[])


# Time to verify if we can compile with xslt
XSLT_CONFIG="xslt-config"
AC_ARG_WITH([xslt-config], 
            [Select a different xslt-config program to detect the configuration options of libxslt.],
            [XSLT_CONFIG=$withval], [AC_CHECK_PROGS([XSLT_CONFIG],[xslt-config],[xslt-config])])

oLIBS="$LIBS"
NEWLIBS="`$XSLT_CONFIG --libs`"
LIBS="$LIBS $NEWLIBS"
oFLAGS="$CFLAGS"
CFLAGS="`$XSLT_CONFIG --cflags` $CFLAGS"
AS_IF([test $USE_XSLT = "true"],
      [AC_CHECK_LIB([xml2], [xmlCreatePushParserCtxt], [
       CFLAGS="$CFLAGS `$XSLT_CONFIG --cflags`"
       AC_DEFINE(HAVE_XSLT)
       USE_XSLT="true"; DUSE_XSLT="-DUSE_XSLT"; OUSE_XSLT=",use_xslt"], [
            LIBS="$oLIBS"
            CFLAGS="$oFLAGS"
            AC_MSG_WARN(I could not find XSLT. I will remove the its support.)
        ])],[ LIBS="$oLIBS"
             CFLAGS="$oFLAGS"
             ])

AS_IF([test $USEINTERFACE = "readline"], [AC_CHECK_LIB([termcap],[tgetent],[LIBS="-ltermcap $LIBS"], [AC_CHECK_LIB([curses],[tgetent],[LIBS="-lcurses $LIBS"], [AC_CHECK_LIB([ncurses],[tgetent],[LIBS="-lncurses $LIBS"], [AC_MSG_FAILURE([The readline interface requires one of the following: termcap, curses, or ncurses.])])])])],[])

AS_IF([test $USEINTERFACE = "readline"], [AC_CHECK_LIB([readline], [readline],
      [LIBS="-lreadline $LIBS"; AC_DEFINE(HAVE_LIBREADLINE)], [AC_MSG_FAILURE([The ncurses interface requires the ncurses library. Either install it, or choose the readline interface. See ./configure --help for more information.])])])

#AS_IF([test $USEPARALLEL = "true"], [AC_CHECK_LIB([$MPILNAME], [MPI_Init],
#      [LIBS="-l$MPILNAME $LIBS"], [AC_MSG_FAILURE([In order to compile the parallel version, an implemention of MPI is required. However, I failed while checking for $MPILNAME.])])])

# Check for malloc, if not there, this is most likely a Mac OS X.
AC_CHECK_HEADER(malloc.h,[],[AC_CHECK_HEADER(/usr/include/malloc/malloc.h,[CFLAGS="-I /usr/include/malloc/ $CFLAGS"], [AC_MSG_FAILURE([Could not find malloc.h])])])

AC_CHECK_PROGS([OCAML], [ocaml])
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
            USE_NATIVEDYNLINK="true"
        fi
	    rm -f str.cmxs*
        ;;
esac

# Final fix to compile in windows
case ${OS_INFORMATION} in
    *-cygwin*)
        LIBS="-link -mwin32 $LIBS"
        ;;
    *)
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

# Checks for library functions.
AC_SEARCH_LIBS([zlibVersion],[z])
AC_SEARCH_LIBS([exp],[m])
AC_PROG_GCC_TRADITIONAL
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_CHECK_FUNCS([floor gettimeofday memset pow sqrt sin exp])

# Time to handle the lapack checks
case ${OS_INFORMATION} in
    *-apple-*)
        AS_IF([test $USE_LIKELIHOOD = "true"],
            [LIBS="-framework veclib $LIBS"; AC_DEFINE(HAVE_LIBLAPACK, 1, "lapack confirmation")],
            [AC_MSG_RESULT([Not enabling Likelihood])])
        ;;
    *)
        AS_IF(  [test $USE_LIKELIHOOD = "false"],
            [AC_MSG_RESULT( [Not enabling Likelihood])],
            [AC_CHECK_LIB([lapack],[main],
                [USE_LIKELIHOOD="true"; DUSE_LIKELIHOOD="-DUSE_LIKELIHOOD"; OUSE_LIKELIHOOD=",use_likelihood"; 
                 AC_MSG_RESULT( [Found lapack and enabling Likelihood]); LIBS="-llapack -lblas -lgfortran $LIBS"],
                [AC_MSG_RESULT( [lapack not found; continuing without Likelihood] )]
        )])
        ;;
esac

CC_VERSION="`$CC --version`"
POY_VERSION_NUMBER="`$MERCURIAL`"

AC_CONFIG_FILES([Makefile _tags myocamlbuild.ml buildNumber.ml compileFlags.ml])
AC_OUTPUT
