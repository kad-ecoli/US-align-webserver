#!/bin/bash
rootdir=`dirname $(readlink -e $0)`
outdir="$rootdir/output"

remove_day=1 # remove after 48 hours

yes|rm $outdir/*pml

cd $outdir

for target in `find US* -mtime +$remove_day`;do
    yes| rm -f $target
done

for prefix in `ls |grep -ohP "US\d{9}"|uniq|sort|uniq`;do
    yes|rm -f `find ${prefix}* -mtime +$remove_day`
done
