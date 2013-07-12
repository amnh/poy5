#!/bin/bash

style=$1
shift
echo "WROTE $style"

cd src/
./configure --enable-interface=ncurses $*

if [$style = "uninstall"];
then
    make uninstall
else
    make poy
    make install
fi
