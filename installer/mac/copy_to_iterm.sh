#!/bin/sh

itermdir="$1"

if test -d $itermdir
then
	for file in `svn list -R`
	do
		if test -f $file
		then
			if test "$2" = "to"
			then
				cp $file $itermdir/$file
			elif test "$2" = "from"
			then
				cp $itermdir/$file $file
			elif test "$2" = "diff"
			then
				diff -u $itermdir/$file $file
			else
				echo 'Expecting second argument to be "to" or "from" '
				exit 1
			fi
		fi
	done
else
	echo "Expecting path to iterm"
	exit 2
fi
