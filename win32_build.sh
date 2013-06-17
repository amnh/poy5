#!/bin/bash

sequential=0
parallel=0
ncurses=0

make_installers=1
update=0

version=${BUILD_NUMBER}

#VARIABLES THAT MIGHT NEED TO BE MODIFIED

    destination="/cygdrive/c/Users/developer/Desktop/poy/${version}"
    concorde="--with-concorde-dir=/usr/i686-w64-minw32/sys-root/mingw/lib/"

    compiler="i686-w64-mingw32-gcc"

    basic_cflags="-I/usr/i686-w64-mingw32/sys-root/mingw/bin -I/usr/i686-pc-mingw32/sys-root/mingw/lib -I/usr/i686-w64-mingw32/sys-root/mingw/include -I/usr/lib/gcc/i686-w64-mingw32/4.5.3"
    basic_lflags="-L/usr/i686-w64-mingw32/sys-root/mingw/bin -L/usr/i686-pc-mingw32/sys-root/mingw/lib -L/usr/i686-w64-mingw32/sys-root/mingw/include -L/usr/lib/gcc/i686-w64-mingw32/4.5.3"

    export FLEXLINKFLAGS="-v -v -L/usr/lib/gcc/i686-mingw32/4.5.3/ -L/cygdrive/c/MPICH2_32/lib -L/usr/i686-w64-mingw32/sys-root/mingw/lib"

    mpi_cflags="-I/cygdrive/c/MPICH2_32/include/ ${basic_cflags}"
    mpi_lflags="-lmpi -L/cygdrive/c/MPICH2_32/lib/ ${basic_lflags}"
    ncurses_cflags=""
    ncurses_lflags="${basic_lflags}"
#
#END

CONFIGURATION_OPTIONS=" --enable-xslt --enable-large-alphabets --enable-long-sequences ${CONFIGURATION_OPTIONS}"

configuration="CC=${compiler} ${concorde} ${CONFIGURATION_OPTIONS} --with-prefix=/usr/i686-w64-mingw32/sys-root/mingw"
while getopts 'uspnm' OPTION; do
    case $OPTION in
        u)
        update=1
        ;;
        s)
        sequential=1
        ;;
        p)
        parallel=1
        ;;
        n)
        ncurses=1
        ;;
        m)
        make_installers=1
        ;;
        ?)
        printf "Usage: \n%s [OPTION]*\n\n-m HOST create installers in the HOST computer (uses ssh for this step)\n-u update from the subversion repository\n-s compile sequential flat \n-p parallel flat version\n-n compile sequential ncurses." $(basename $0)
        exit 2
        ;;
    esac
done

function compile_executable {
    make clean # delete crap because ocamlbuild on windows is STUPID
    if ! ocamlbuild poy.native -lflags "-nodynlink -cclib -link -cclib -static" -ocamlopt "ocamlopt.opt -verbose"; then
        echo "I could not make poy!!! ..."
        exit 1
    fi
}

if [ $update -eq 1 ]; then
    if ! hg pull; then
        echo "Repository pull failed ... sorry pal!"
        exit 1
    fi
    if ! hg update "${TAG_NUMBER}"; then
        echo "Repository update failed ... sorry pal!"
        exit 1
    fi
fi

rm -rf ncurses_poy.exe
rm -rf par_poy.exe
rm -rf seq_poy.exe

rm -rf $destination
mkdir $destination

cd src
if [ $ncurses -eq 1 ]; then
    # We first compile the regular ncurses interface
    if ! ./configure $configuration --enable-xslt --enable-interface=pdcurses CFLAGS="$ncurses_cflags" LIBS="$ncurses_lflags"; then
        echo "Configuration failed!"
    else
        compile_executable
        if [ -e _build/poy.native ]; then
            cp _build/poy.native ../ncurses_poy.exe
        fi
    fi
    if [ $make_installers -eq 1 ]; then
       if ! cp -f ./_build/poy.native ${destination}/ncurses_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   else
       if ! cp -f ./_build/poy.native ncurses_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   fi
fi

if [ $sequential -eq 1 ]; then
    # Now we compile the html interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html CFLAGS="$basic_cflags" LIBS="$basic_lflags"; then
        echo "Configuration failed!"
    else
        compile_executable
        if [ -e _build/poy.native ]; then
            cp _build/poy.native ../seq_poy.exe
        fi
    fi
   if [ $make_installers -eq 1 ]; then
       if ! cp -f ./_build/poy.native ${destination}/seq_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   else
       if ! cp -f ./_build/poy.native seq_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   fi
fi

if [ $parallel -eq 1 ]; then
    # Now we compile the parallel interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html --enable-mpi CFLAGS="$mpi_cflags" LIBS="$mpi_lflags"; then
        echo "Configuration failed!"
    else
        compile_executable
        if [ -e _build/poy.native ]; then
            cp _build/poy.native ../par_poy.exe
        fi
    fi
    if [ $make_installers -eq 1 ]; then
       if ! cp -f ./_build/poy.native ${destination}/par_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   else
       if ! cp -f ./_build/poy.native par_poy.exe; then
           echo "I could not replace the executable in the distribution"
           exit 1
       fi
   fi
fi

#if [ $make_installers -eq 1 ]; then
#    rm -f /cygdrive/c/POY_Installer.exe
#    ./create_installers.bat
#    if ! scp /cygdrive/c/POY_Installer.exe ${MACHOST}:poy_build/poy/${version}/windows/; then
#        echo "I could not copy the resulting executable in $MACHOST!"
#        exit 1
#    fi
#fi
