#!/bin/bash


prc='runrsQA.C'
if [ $# -gt 0 ] ; then
#  ddd  aliroot  $1
    prc=$1
fi

aliroot -b -n -q $prc >& proc.log



