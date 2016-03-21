#!/bin/bash
##################################################
validateout=`dirname $0`
validatetime=`date`
validated="0";
error=0
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout
echo "* Time:    $validatetime " >> stdout
echo "* Dir:     $validateout" >> stdout
echo "* Workdir: $validateworkdir" >> stdout
echo "* ----------------------------------------------------*" >> stdout
ls -la ./ >> stdout
echo "* ----------------------------------------------------*" >> stdout

##################################################


if [ ! -f stderr ] ; then
   error=1
   echo "* ########## Job not validated - no stderr  ###"  >> stdout
   echo "Error = $error"  >> stdout
fi
segViol=`grep -Ei "Segmentation violation" stderr`
segFault=`grep -Ei "Segmentation fault" stderr`
glibcErr=`grep -Ei "*** glibc detected ***" stderr`

if [ "$segViol" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. violation  ###"  >> stdout
   echo "$segViol"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$segFault" != "" ] ; then
   error=1
   echo "* ########## Job not validated - Segment. fault  ###"  >> stdout
   echo "$segFault"  >> stdout
   echo "Error = $error"  >> stdout
fi
if [ "$glibcErr" != "" ] ; then
   error=1
   echo "* ########## Job not validated - *** glibc detected ***  ###"  >> stdout
   echo "$glibcErr"  >> stdout
   echo "Error = $error"  >> stdout
fi

if [ $error = 0 ] ; then
   echo "* ----------------   Job Validated  ------------------*" >> stdout
fi
echo "* ----------------------------------------------------*" >> stdout
echo "*******************************************************" >> stdout
cd -
exit $error
