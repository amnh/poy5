#!/bin/sh
mkdir binaries
rm -f binaries/*
rm -f mpoy.dmg
cp mpoy-* binaries/
hdiutil create -srcfolder binaries mpoy.dmg
rm -fr binaries
