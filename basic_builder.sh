#!/bin/bash

cd src/
./configure --enable-interface=ncurses $@
make poy
make install
