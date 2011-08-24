#!/bin/bash
#
# This script is used for building the system remotely in our build system. It
# does contain a basic outline of compilation, but requires a release number
# passed to it (although will work fine without it), will statically link the
# application, and makes assumptions of where the static libraries of LAPACK and
# BLAS are located. Result from command is placed in {this_directory}/linux/.
#
# Arguments 
#   1 - release number 

export FLAGS="-static -static-libgcc -O3 -msse3"

LINUX_DIRECTORY=linux
rm -rf linux
mkdir  linux


# Build the version to be paired with the GUI
if ! ./configure --with-version-number=$1 --enable-interface=html CFLAGS="${FLAGS}" ; then
    echo "Failure in html interface configuration"
    exit 1
fi

if ! make clean; then
    echo "Could not clean!"
    exit 1
fi
if ! make poy.native; then
    echo "Failure in make step"
    exit 1
fi
cp ./src/_build/poy.native ./$LINUX_DIRECTORY/seq_poy.command

# Now we make the ncurses interface
if ! ./configure --with-version-number=$1 --enable-interface=ncurses CFLAGS="${FLAGS}" ; then
    echo "Failure in ncurses interface configuration"
    exit 1
fi
if ! make clean; then
    echo "Could not clean!"
    exit 1
fi
if ! make poy.native; then
    echo "Failure in make step"
    exit 1
fi

cp ./src/poy.native ./$LINUX_DIRECTORY/ncurses_poy
#package ncurses around xterm call
cat > ./${LINUX_DIRECTORY}/ncurses_poy.command <<EOF
#!/bin/bash
xterm -e ../Resources/ncurses_poy
EOF
chmod a+x ./${LINUX_DIRECTORY}/ncurses_poy.command
