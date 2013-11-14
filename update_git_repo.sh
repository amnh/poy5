#!/usr/bin/env bash
#
# We would like to mirror on github, but we have files that are large in our
# history that cannot be represented in github. This small script will prune the
# test/ directory in our distribution thus allowing mirroring. This will remove
# all the changes at the amnh/poy repo on github so contributions need to be
# mirrored on mercurial first.
#
# This script will take awhile --approx 30minutes to complete.
#
# usage : bash update_git_repo.sh source_repo destination_repo github_repo
#
SRC=${1}
DEST=${2}

# clean source
cd ${SRC}
hg clean
# clean destination; repo must be !!bare!!
rm -r ${DEST}
git init --bare ${DEST}
# push to destination
hg push ${DEST}
#prune destination
cd ${DEST}
git filter-branch --force --index-filter 'git rm -rf --cached --ignore-unmatch test' --prune-empty --tag-name-filter cat -- --all
#push to github
git remote add origin ssh://git@github.com/amnh/POY
git push --force -u origin master
