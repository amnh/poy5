#!/bin/bash
# USAGE: ./bisect.sh <poy script> <lowerbound-good;or 773> <upperbound-bad; or current>

POYSCRIPT=$1
CHANGESET1=$2
CHANGESET2=$3
LOG="bisect.log"

TESTING=0   # set to 1 to skip modifcation check for testing purposes

if [ -z "$CHAGESET1" ];
then 
    CHANGESET1=773
fi

#ensure that the repo has no modifications
if hg status -m | grep -q "[MA]"
then
    if [ $TESTING -eq 0 ] 
    then
        echo "Please push/revert changes in current revision."
        exit 0
    fi
fi

#clear log file
echo > $LOG
# reset the bisect functionality
hg bisect --reset
# set the lower bound; good
hg bisect --good $CHANGESET1
#set the upper bound; bad,
hg bisect --bad $CHANGESET2

while :
do
    CURRENT_GLOBAL=`hg identify | sed -e 's/ .*//g'`
    CURRENT_LOCAL=`hg identify -n | sed -e 's/ .*//g'`
    # clean build
   ocamlbuild -clean
    # configure; pipe to null; we don't need to see this crap
   ./configure --enable-interface=readline --enable-likelihood > /dev/null
    # turn on debugging options in certain files. also fpic in early builds.
   sed -i -e 's/debug_pass_errors \= false/debug_pass_errors \= true/' poy.ml
   sed -i -e 's/let cflags = "\(.*\?\)"/let cflags = "-fPIC \1"/' myocamlbuild.ml
    # compile
   ocamlbuild poy.native
    # run
   ./poy.native $POYSCRIPT
    # revert previous file; myocamlbuild is generated and can be ignored
   hg revert poy.ml
    # ask user if this changeset is correct, to skip or quit, and bisect based on results
   read -p "Did the script work? (y/n/s|q):  "
   [ "$REPLY" != "y" ] || TYPE='-g'
   [ "$REPLY" != "n" ] || TYPE='-b'
   [ "$REPLY" != "s" ] || TYPE='-s'
   [ "$REPLY" != "q" ] || exit 0    # bisect can be proceeded manually if desired

   # bisect; and test the output in case we have finished
   if 
       hg bisect $TYPE | grep -q 'Testing changeset'
   then
       echo $CURRENT_LOCAL $CURRENT_GLOBAL $TYPE >> $LOG
   else
       # report final information to user
       hg summary
       hg summary >> $LOG
       exit 0
   fi
done
