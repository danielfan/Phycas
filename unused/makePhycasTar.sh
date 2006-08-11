#!/bin/sh
if [ $# -lt 2 ]; then
	echo "2 arguments are required:"
	echo "   $0 <path_to_top_dir> <tarball-prefix> "
	exit
fi
if [ ! -x $1 ]; then
	echo "The top level directory ($1) was not found"
	exit 1
fi
tarParent="$2"
if [ -x $tarParent ] || [ -f $tarParent ]; then
	echo "The  file or directory $tarParent is in the way."
	exit 2
fi
`mkdir $tarParent` || exit
for subdir in cipres_services phycas ncl ; do
	`cd $tarParent; cvs checkout phycasdev/$subdir`
	`mv $tarParent/phycasdev/$subdir $tarParent/`
done
`rm -r $tarParent/phycasdev`
for symLinkToCopy in missing depcomp install-sh config.sub config.guess ; do
	`cp $1/$symLinkToCopy $tarParent/`
done
for symLinkToCopy in COPYING INSTALL ; do
	`cp $1/cipres_services/$symLinkToCopy $tarParent/cipres_services`
done
#`cd $1; set PYTHONPATH=(./$1/python); export PYTHONPATH;  python freshmake.py --code;`
`cp -r $1/phycas/command $tarParent/phycas/` || exit
`cd $tarParent/cipres_services; sh reconf`
`tar czf $tarParent.tgz $tarParent`
