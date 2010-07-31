


if ! test -f setuptools-0.6c11-py2.7.egg
then
    wget http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg#md5=fe1f997bc722265116870bc7919059ea
fi
sh setuptools-0.6c11-py2.7.egg --prefix=$SELF_CONTAINED_PREFIX



if ! test -d ipython-0.10
then
    wget http://ipython.scipy.org/dist/0.10/ipython-0.10.tar.gz
    tar xfvz ipython-0.10.tar.gz
fi
cd ipython-0.10 || exit
python setup.py install --prefix=$SELF_CONTAINED_PREFIX || exit
cd .. || exit

cp  trunk_phycas/installer/mac/iterm_with_python/iTerm4PhycasexecFromGUI.sh PhycasGUI.app/Contents/MacOS/iTerm4PhycasexecFromGUI.sh || exit
