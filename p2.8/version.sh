#!/bin/sh
#
# version.sh: writes revision and date to version.h
#
# $Date$
# $Revision$
# $Id$

rm -f version.h

date=`date +"%Y-%m-%d"`

if [ -d "../.git" ];
then
	revision=`git describe --always`;
else
	revision="NotGitRepo";
fi
current_folder=`pwd`

compilation_date=`date`

echo \#define DATE \"$date\" > version.h
echo \#define REVISION \"$revision\" >> version.h
echo \#define COMPILE_DATE \"$compilation_date\" >> version.h
echo \#define LOCATION \"$current_folder\" >> version.h

echo $revision
