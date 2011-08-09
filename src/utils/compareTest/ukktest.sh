#this script create random int array of random length(1-100). as input of two
#functions -- here they are space saving version of ukkonen and full alignment
#version, compare the result -- here they are the cost of alignment.
#if two results are the same, print them to screen, this test pass, 
# continue with next test. other wise , test fail.

#!/bin/bash
while true; do
ocamlopt -o randomarr.native randomarr.ml
ocamlopt -o ukkss.native ukk_space_save.ml
#ocamlopt -o ukkf.native ukkfloat.ml
ocamlopt -o algnfull.native algnfull.ml
./randomarr.native 
#./ukkf.native > ukkf.out2
./algnfull.native > ukkf.out2
./ukkss.native > ukkss.out2
read var1 < ukkss.out2
echo "var1=$var1"
read var2 < ukkf.out2
echo "var2=$var2"
if [ $var1 -eq $var2 ]; then
    echo "Pass"
else
    echo "Fail"
    exit 1
fi
rm ukkf.out2 ukkss.out2
done

#out1=./ukkss.native
#out2=./ukkf.native
#$1 2>&1 > $out1
#$2 2>&1 > $out2
#diff -q >/dev/null $out1 $out2
#if [ $? != 0 ]; then
#    echo "FAIL"
#    echo
#    echo "Details:"
#    diff $out1 $out2
#    echo
#    exit 1
#else
#    echo "PASS"
#fi
#rm $out1 $out2
#done
#echo "diff of result:"
#diff <(./ukkss.native) <(./ukkf.native)
#echo "common of result:"
#comm <(./ukkss.native) <(./ukkf.native)
#if [ $S1 = $S2 ];
#then 
#    echo "S1('$S1') is equal to S2('$S2')"
#fi
