#!/bin/bash
# This script takes care of generating binaries for all mac plataforms.

# A function that takes four arguments: the SDK to use, the target architecture,
# the configuration options, and finally the target binary where poy will be
# placed.
UNIVERSAL_DIRECTORY=universal

function generate_binary {
    echo "Configurating for $1 - $2 - $3 for target $4" >> distro.log
    echo "Configure for $1 - $2 - $3 for target $4" 
    if ! PATH=$5:$PATH CFLAGS="-O3 -arch $2 -isysroot $1 $6" CC="$7" ./configure $3 >> distro.log; then
        echo "Failed in the configuration step for $1 $2 $3 target $4"
        exit 1
    fi
    echo "Make Clean"
    if ! PATH=$5:$PATH make clean >> distro.log; then
        echo "Failed in the clean step for $1 $2 $3 target $4"
        exit 1
    fi
    if [ ! $8 = "" ]; then
        echo "Make OcamlMPI"
        PATH=$5:$PATH make ocamlmpi >> distro.log
        echo "Make Depend"
        if ! PATH=$5:$PATH make depend >> distro.log; then
            echo "Failed in the depend step for $1 $2 $3 target $4"
            exit 1
        fi
        echo "Make $8"
        cd src
    fi
    if ! PATH=$5:$PATH make $8 >> distro.log; then
        echo "Failed in the make step for $1 $2 $3 target $4"
        cd ../
        exit 1
    fi
    if [ ! $8 = "" ]; then
        cd ../
        mv -f src/$8 $4
    else
        mv -f src/poy $4
    fi
}

OCAMLROOT=/opt/ocaml
OCAMLVERSION=3.10.0
OCAML_PATH=/opt/ocaml/${OCAMLVERSION}

# Generating for Panther, we don't produce parallel version.
function panther_distribution {
    rm -f ./panther/seq_poy.command
    rm -f ./panther/ncurses_poy.command
    if ! generate_binary /Developer/SDKs/MacOSX10.3.9.sdk ppc \
        "--enable-xslt --enable-interface=html" ./panther/seq_poy.command \
        ${OCAML_PATH}/panther/bin "" gcc ""; then
        exit 1
    fi

    if ! generate_binary /Developer/SDKs/MacOSX10.3.9.sdk ppc \
        "--enable-xslt --enable-interface=ncurses" ./panther/ncurses_poy.command \
        ${OCAML_PATH}/panther/bin  "" gcc ""; then
        exit 1
    fi
}

# Generating for universal binaries under Tiger.

# Generating for PowerPC
function generate_ppc {
    generate_binary /Developer/SDKs/MacOSX10.4u.sdk ppc "$1" \
    ./$UNIVERSAL_DIRECTORY/$3 ${OCAML_PATH}/tiger/ppc/bin "" "$2" $4
}

function generate_intel {
    generate_binary /Developer/SDKs/MacOSX10.4u.sdk i386 "$1" \
    ./$UNIVERSAL_DIRECTORY/$3 ${OCAML_PATH}/tiger/intel/bin "" "$2" $4
}

function make_universals {
    lipo -output ./$UNIVERSAL_DIRECTORY/$1 -create ./$UNIVERSAL_DIRECTORY/$2 \
    ./$UNIVERSAL_DIRECTORY/$3
}

function tiger_distribution {
    if ! generate_ppc "--enable-xslt --enable-interface=html" gcc seq_poy_ppc ""; then
        exit 1
    fi

    if ! generate_ppc "--enable-xslt --enable-interface=ncurses" gcc ncurses_poy_ppc ""; then
        exit 1
    fi

    if ! generate_ppc "--enable-xslt --enable-interface=html --enable-mpi" \
        /usr/local/poy4/mpich2-1.0.5p2/gforker/arch/bin/mpicc par_poy_pcc ""; then 
        exit 1
    fi

    if ! generate_intel "--enable-xslt --enable-interface=html" gcc seq_poy_intel ""; then
        exit 1
    fi

    if ! generate_intel "--enable-xslt --enable-interface=ncurses" gcc ncurses_poy_intel ""; then
        exit 1
    fi

    if ! generate_intel "--enable-xslt --enable-interface=html --enable-mpi" \
        /usr/local/poy4/mpich2-1.0.5p2/gforker/arch/bin/mpicc_icc \
        par_poy_intel ""; then
        exit 1
    fi

    # Finally, generate the universal binaries
    if ! make_universals seq_poy.command seq_poy_ppc seq_poy_intel; then
        echo "Failed in lipo sequential for universal binaries"
        exit 1
    fi

    if ! make_universals par_poy.command par_poy_pcc par_poy_intel; then
        echo "Failed in lipo parallel for universal binaries"
        exit 1
    fi

    if ! make_universals ncurses_poy.command ncurses_poy_ppc ncurses_poy_intel; then
        echo "Failed in lipo ncurses for universal binaries"
        exit 1
    fi
}

function create_generic {
    # We want to be able to create universal binaries for any target basically.
    if ! generate_intel "$1" /opt/intel/cc/9.1.037/bin/icc intel_$2 $2; then
        echo "Failed generating the intel part of the target $2";
        exit 1
    fi
    if ! generate_ppc "$1" gcc ppc_$2 $2; then
        echo "Falied generating the ppc part of the target $2";
        exit 1;
    fi
    if ! make_universals $2 intel_$2 ppc_$2; then
        echo "Falied in lipo for $2"
        exit 1;
    fi
}

target=
configuration=
while getopts 't:c:h:help' OPTION; do
    case $OPTION in
        t) 
        target="$OPTARG"
        ;;
        c) 
        configuration="$OPTARG"
        ;;
        ?)
        printf "Usage: \n%s -t [panther | tiger | target] -c \"configuration_options\"\n" $(basename $0)
        exit 2
    esac
done

echo "Building target $target"

if [ "panther" == $target ]; then
   panther_distribution
elif [ "tiger" == $target ]; then
   tiger_distribution
else
    rm -f $UNIVERSAL_DIRECTORY/$target
    create_generic "$configuration" $target
fi
