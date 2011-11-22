#!/bin/bash

sequential=0
parallel=0
ncurses=0
make_installers=0
update=0
version=${BUILD_VERSION}

#VARIABLES THAT MIGHT NEED TO BE MODIFIED
#
#   # FOR OUR MYSTERIOUS WINDOWS XP MACHINE
#   destination="/cygdrive/c/poy_distribution/"
#   concorde="--with-concorde-dir=/home/andres/concorde"
#   basic_cflags="-msse3 -O3 -I/home/andres/zlib/include" 
#   basic_lflags="-L/home/andres/zlib/lib"
#   mpi_cflags="-I/cygdrive/c/mpich2/include ${basic_cflags}"
#   mpi_lflags="-L/cygdrive/c/mpich2/lib -lmpi ${basic_lflags}"
#   ncurses_cflags="-I/home/andres/PDCurses-3.3/ ${basic_cflags}"
#   ncurses_lflags="-L/home/andres/PDCurses-3.3/win32/ ${basic_lflags}"
#
    # FOR OUR SLIGHTLY LESS MYSTERIOUS WINDOWS 7 MACHINE(S)
    destination="/cygdrive/c/poy_builds/${BUILD_VERSION}"
    concorde="--with-concorde-dir=/home/Developer/programs/concorde/"
    basic_cflags="-msse3 -O3 -I/usr/i686-pc-mingw32/xslt_xml2/include/ -L/usr/i686-pc-mingw32/xslt_xml2/lib/"
    basic_lflags="${LFLAGS}"
    mpi_cflags="-I/cygdrive/c/MPICH2/include/ ${basic_cflags}"
    mpi_lflags="-L/cygdrive/c/MPICH2/lib/ -lmpi ${basic_lflags}"
    ncurses_cflags="-I/home/Developer/programs/PDCurses-3.4/ -L/home/Developer/programs/PDCurses-3.4/mingw-build/ ${basic_cflags}"
    ncurses_lflags="${basic_lflags}"
#
#END

configuration="${concorde} ${CONFIGURATION_OPTIONS} --with-prefix=/usr/i686-pc-mingw32/sys-root/mingw"
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
    bash clean.sh # delete crap because ocamlbuild on windows is STUPID
    if ! ocamlbuild poy.native; then
        echo "I could not make poy!!! ..."
        exit 1
    fi
}

if [ $update -eq 1 ]; then
    if ! hg pull; then
        echo "Repository pull failed ... sorry pal!"
        exit 1
    fi
    if ! hg update ${TAG_NUMBER}; then
        echo "Repository update failed ... sorry pal!"
        exit 1
    fi
fi

rm -rf ${destination}

mkdir ${destination}

cd src

if [ $ncurses -eq 1 ]; then
    # We first compile the regular ncurses interface
    if ! ./configure $configuration --enable-xslt --enable-interface=ncurses CFLAGS="$ncurses_cflags" LFLAGS="$ncurses_lflags"; then
        echo "Configuration failed!"
        exit 1
    fi
    compile_executable
    cp _build/poy.native ../ncurses_poy.exe
#   if [$make_installers -eq 1]; then
#       if ! cp -f ./_build/poy.native ${destination}/bin/ncurses_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   else
#       if ! cp -f ./_build/poy.native ncurses_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   fi
fi

if [ $sequential -eq 1 ]; then
    # Now we compile the html interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html CFLAGS="$basic_cflags" LFLAGS="$basic_lflags"; then
        echo "Configuration failed!"
        exit 1
    fi
    compile_executable
    cp _build/poy.native ../seq_poy.exe
#   if [ $make_installers -eq 1 ]; then
#       if ! cp -f ./_build/poy.native ${destination}/bin/seq_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   else
#       if ! cp -f ./_build/poy.native seq_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   fi
fi

if [ $parallel -eq 1 ]; then
    # Now we compile the parallel interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html --enable-mpi CFLAGS="$mpi_cflags" LFLAGS="$mpi_lflags"; then
        echo "Configuration failed!"
        exit 1
    fi
    compile_executable
    cp _build/poy.native ../par_poy.exe
#   if [$make_installers -eq 1]; then
#       if ! cp -f ./_build/poy.native ${destination}/bin/par_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   else
#       if ! cp -f ./_build/poy.native par_poy.exe; then
#           echo "I could not replace the executable in the distribution"
#           exit 1
#       fi
#   fi
fi

#if [ $make_installers -eq 1 ]; then
#    rm -f /cygdrive/c/POY_Installer.exe
#    ./create_installers.bat
#    if ! scp /cygdrive/c/POY_Installer.exe ${MACHOST}:poy_build/poy/${version}/windows/; then
#        echo "I could not copy the resulting executable in $MACHOST!"
#        exit 1
#    fi
#fi
