#!/bin/bash
# The environment MACHOST holds the URL of the macintosh host running the
# virtual machine. The only argument contains the desired extra flags for the
# configuration script (for example long sequence support). 
export PATH=/home/andres/minglibs/bin:/cygdrive/c/ocamlmgw/3_11_0/bin:$PATH
source $HOME/.keychain/${HOSTNAME}-sh

sequential=0
parallel=0
ncurses=0
configuration=${CONFIGURE_OPTIONS}
MACHOST=""
make_installers=0
update=0
version=${BUILD_VERSION}

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
        MACHOST="samson"
        ;;
        ?)
        printf "Usage: \n%s [OPTION]*\n\n-m HOST create installers in the HOST computer (uses ssh for this step)\n-u update from the subversion repository\n-s compile sequential flat \n-p parallel flat version\n-n compile sequential ncurses." $(basename $0)
        exit 2
        ;;
    esac
done

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

function compile_executable {

cd ./src

if ! make clean; then
    echo "I could not clean up the distribution ..."
    exit 1
fi

if ! make poy.native; then
    echo "I could not make poy!!! ..."
    exit 1
fi

cd ../
}

if [ $ncurses -eq 1 ]; then
    # We first compile the regular ncurses interface
    if ! ./configure $configuration --enable-xslt --enable-interface=ncurses CFLAGS="-O3 -msse3 -L/home/andres/zlib/lib -L/home/andres/PDCurses-3.3/win32/ -I /home/andres/zlib/include -I /home/andres/PDCurses-3.3/" LIBS="-L/home/andres/zlib/lib -L/home/andres/PDCurses-3.3/win32/" ; then
        echo "Configuration failed!"
        exit 1
    fi
    compile_executable
    if ! cp -f ./src/poy.native /cygdrive/c/poy_distribution/bin/ncurses_poy.exe; then
        echo "I could not replace the poy executable in the distribution"
        exit 1
    fi
fi

if [ $sequential -eq 1 ]; then
    # Now we compile the html interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html CFLAGS="-msse3 -O3 -L/home/andres/zlib/lib  -I /home/andres/zlib/include" LIBS="-L/home/andres/zlib/lib -L/home/andres/PDCurses-3.3/win32/ "; then
        echo "Configuration failed!"
        exit 1
    fi
    compile_executable
    if ! cp -f ./src/poy.native /cygdrive/c/poy_distribution/bin/seq_poy.exe; then
        echo "I could not replace the executable in the distribution"
        exit 1
    fi
fi

if [ $parallel -eq 1 ]; then
    # Now we compile the parallel interface
    if ! ./configure $configuration --enable-xslt --enable-interface=html --enable-mpi CFLAGS="-msse3 -O3 -L/home/andres/zlib/lib -L/cygdrive/c/mpich2/lib -I /home/andres/zlib/include -I /cygdrive/c/mpich2/include"  LIBS="-L/home/andres/zlib/lib -L/home/andres/PDCurses-3.3/win32/ -L/cygdrive/c/mpich2/lib -lmpi"; then
        echo "Configuration failed!"
        exit 1
    fi
    make clean
    make
    if ! cp -f ./src/poy.native /cygdrive/c/poy_distribution/bin/par_poy.exe; then
        echo "I could not replace the executable in the distribution"
        exit 1
    fi
fi

if [ $make_installers -eq 1 ]; then
    rm -f /cygdrive/c/POY_Installer.exe
    ./create_installers.bat
    if ! scp /cygdrive/c/POY_Installer.exe ${MACHOST}:poy_distro/source_code/binaries/windows/; then
        echo "I could not copy the resulting executable in $MACHOST!"
        exit 1
    fi
fi
