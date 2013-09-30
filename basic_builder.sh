#!/bin/bash

style=$1
shift

cd src/
./configure $*

if [ $style = "uninstall" ];
then
    make uninstall
else
    make poy
    make install
fi
