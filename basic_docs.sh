#!/bin/bash

cd src/
./configure --enable-interface=ncurses $@
make docs
make install
