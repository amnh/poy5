#!/bin/bash
# USAGE: ./bisect.sh <poy script> <lowerbound-good;or 773> <upperbound-bad; or current>

POYSCRIPT=$1
CHANGESET1=$2
CHANGESET2=$3

if [$CHAGESET1 = ""] 
then CHANGESET1=773 
fi

#ensure that the repo has no modifications
if 
    hg status -m 
then
    echo "Please push/revert changes in current revision."
    exit 0
fi

# reset the bisect functionality
hg bisect --reset
# set the lower bound; good
hg bisect --good $CHANGESET1
#set the upper bound; bad, 
hg bisect --bad $CHANGESET2

while :
do
    # clean build
   ocamlbuild -clean
    # configure
   ./configure --enable-interface=readline --enable-likelihood
    # turn on debugging option
   sed -i -e 's/debug_pass_errors \= false/debug_pass_errors \= true/' poy.ml
   sed -i -e 's/let cflags = "\(.*\?\)"/let cflags = "-fPIC \1"/' myocamlbuild.ml
    # compile
   ocamlbuild poy.native
    # run
   ./poy.native $POYSCRIPT
    # revert previous file
   hg revert poy.ml
    # ask user if this changeset is correct, and bisect based on results
   read -p "Did the script work? (y/n/s):  "
   [ "$REPLY" != "y" ] || TYPE='-g' 
   [ "$REPLY" != "n" ] || TYPE='-b'
   [ "$REPLY" != "s" ] || TYPE='-s'
    
   # bisect; and test the output in case we have finished
   if 
       hg bisect $TYPE | grep -q 'Testing changeset'
   then 
       echo
   else
       # report final information to user
       hg summary
       exit 0
   fi

done
