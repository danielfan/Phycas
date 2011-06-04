#!/usr/bin/env python
'''A script to automate the building of phycas.

Expects the current directory structure to be:
    - Phycas (or $PHYCAS_ROOT)
    - boost_1_4*** (or $BOOST_ROOT)
    - ncl (or $NCL_INSTALL_DIR)
'''

import sys
import os
import subprocess
import urllib
import tarfile
import platform

def untar(fn):
    try:
        tar = tarfile.open(fn)
    except:
        die("The tar archive "+  fn + " could not be opened.")
    try:
        tar.extractall()
        tar.close()
    except:
        die("The tar archive "+  fn + " could not be unpacked. It may be corrupted.\n\nTry removing the file and running this script again. \nYou may need to remove the directory created by the partial unpacking of this archive, as well.")
    return True

################################################################################
# Globals
################################################################################

SCRIPT_NAME = os.path.split(sys.argv[0])[-1]
IS_MAC = platform.uname()[0].lower().startswith('darwin')
try:
    NUM_PROC = str(int(os.environ.get('NUM_CORES_FOR_BUILD')))
except:
    NUM_PROC = '1'
################################################################################
# UI helpers
################################################################################
def debug(*valist):
    msg = ' '.join([str(i) for i in valist])
    sys.stderr.write(SCRIPT_NAME + ': ' + msg + '\n')

status = debug
warn = status

def die(*valist):
    warn(*valist)
    sys.exit(1)
    
def positive_response(prompt):
    resp = raw_input(prompt + ' (y/n) ')
    return resp.lower().startswith('y')

def choice_prompt(prompt, choices, num_tries=1):
    n = 0
    while True:
        resp = raw_input(prompt + ' Your choices:\n    ' + '\n    '.join(choices) + '\n')
        n += 1
        if n < num_tries and resp.lower() not in [i.lower() for i in choices]:
            warn(resp, 'is not a valid choice.\n')
        else:
            return resp

################################################################################
# Helpers
################################################################################
def download_status_hook(a,b,c): 
    "From Xuan's response at http://stackoverflow.com/questions/974741/wget-vs-urlretrieve-of-python"
    sys.stderr.write("% 3.1f%% of %d bytes\r" % (min(100, float(a * b) / c * 100), c))
    sys.stderr.flush()

def do_http_download(url, file):
    status('Downloading from ', url)
    urllib.urlretrieve(url, file, download_status_hook)
    sys.stderr.write('\n')
    return True

def do_system(invoc, dir=None):
    if dir:
        orig = os.path.abspath(os.curdir)
    try:
        if dir:
            os.chdir(dir)
        debug('Executing command:\n   ', ' '.join([repr(i) for i in invoc]), '\nin', os.path.abspath(os.curdir))
        try:
            retcode = subprocess.call(invoc, shell=False)
            if retcode == 0:
                return True
            warn(invoc[0] + "died with error code", retcode)
        except OSError, e:
            warn(invoc[0] + " execution failed:", str(e))
        return False
    finally:
        if dir:
            os.chdir(orig)
    
def do_svn_checkout(url, filepath):
    return do_system(["svn", "checkout", url, filepath])

def do_git_checkout(url, filepath):
    return do_system(["git", "clone", url, filepath])

wrong_py_version_msg = 'You will need to install Python version 2.4.*, 2.5.*, 2.6.*, or 2.7.* to build Phycas.\nRerun this script after installing the appropriate python'
try:
    import platform
    py_major, py_minor, py_patch = [int(i) for i in platform.python_version().split('.')]
except:
    die("Calls to python's platform module failed. ", wrong_py_version_msg)
debug("Python version", py_major, ".", py_minor, ".", py_patch, "detected")
if py_major != 2 or py_minor not in [4, 5, 6, 7]:
    die(wrong_py_version_msg)





################################################################################
# Dependencies
################################################################################

class CodeResource(object):
    def __init__(self,
                 name,
                 protocol,
                 url,
                 compressed_name=None,
                 uncompress_fn=None, 
                 uncompressed_name=None,
                 version=None
                 ):
        self.name = name
        self.protocol = protocol.lower()
        self.url = url
        self.compressed_name = compressed_name
        self.uncompress_fn = uncompress_fn
        self.uncompressed_name = uncompressed_name
        self.version = version

    def download_and_unpack(self):
        url = self.url
        file = self.compressed_name
        unpack_cmd = self.uncompress_fn
        final_dir = self.uncompressed_name
        protocol = self.protocol
        
        if os.path.exists(final_dir):
            warn(final_dir, ' already exists and is in the way of the download operation')
            return False
        if file and os.path.exists(file):
            status(file, ' found and will be used instead of download.')
        else:
            if protocol == 'http':
                if do_http_download(url, file):
                    status(file, 'downloaded.')
                else:
                    warn('Download from ', url, 'failed')
                    return False
            elif protocol == 'svn':
                if do_svn_checkout(url, final_dir):
                    return True
            elif protocol == 'git':
                if do_git_checkout(url, final_dir):
                    return True
            else:
                die('Unknown protocol', protocol)
                
        if unpack_cmd:
            if os.path.exists(file):
                status('Unpacking', file, '...')
                if not unpack_cmd(file):
                    die('Could not unpack', file)
        return os.path.exists(final_dir)

class Dependency(object):
    def __init__(self, name, env_var, search_dir_prefix, resource_list):
        self.name = name 
        self.env_var = env_var
        self.search_dir_prefix = search_dir_prefix
        self.resource_list = resource_list
        self.full_path = None
    def acquire(self):
        '''Get the dependency using the following cascade:
            1. env var
            2. subdir that starts with search_dir_prefix
            3. user prompted path
            4. download from resource list
        '''
        env_value = os.environ.get(self.env_var)
        approved_download = False
        if env_value:
            env_value = os.path.abspath(env_value)
            debug(self.env_var, "is set to", env_value)
            if not (os.path.exists(env_value)):
                warn("The path", env_value, "does not exist.")
                env_value = ''
            elif not os.path.isdir(env_value):
                warn("The path", env_value, "is not a directory.")
                env_value = ''
            elif not positive_response("Do you want to use this as the basis for your " + self.env_var + " variable during your Phycas build?"):
                warn("Ignoring", self.env_var)
                env_value = ''
        if not env_value:
            subdirs = [i for i in os.listdir(os.curdir) if os.path.isdir(i)]
            poss_env_value_list = [d for d in subdirs if d.lower().startswith(self.search_dir_prefix)]
            if len(poss_env_value_list) == 1:
                d = poss_env_value_list[0]
                if positive_response(d + ' found. Would you like to use this as your ' + self.env_var+ '?'):
                    env_value = os.path.abspath(d)
            else:
                if len(poss_env_value_list) > 0:
                    d = choice_prompt('Multiple ' + self.env_var + ' directories were found. Enter the name of the one that you would like to use.',
                                     poss_env_value_list,
                                     num_tries=1)
                else:
                    d = raw_input("Enter 'y' to download " + self.name + ", or enter the path to the directory that you would like to use as " + self.env_var + ": ")
                    if d.lower() == 'y':
                        d = ''
                        approved_download = True
                if d:
                    if not os.path.isdir(d):
                        die(d, 'is not a valid directory')
                    env_value = os.path.abspath(d)
            
        if env_value:
            self.full_path = os.path.abspath(env_value)
            return True
    
        if approved_download or positive_response("Would you like to download " + self.name + "?"):
            for resource in self.resource_list:
                if resource.download_and_unpack():
                    self.full_path = os.path.abspath(resource.uncompressed_name)
                    return True
                warn("Attempt to obtain " + self.name + " via " + resource.protocol + " failed. Moving on to next method...")
            die("Could not obtain " + self.name + " libraries. Exiting...")
        die("Could not locate " + self.name + " libraries. Exiting...")

def write_sourceable_file(o, env_dict):
    o.write('''#!/bin/sh\n''')
    for opt_key in ['CXXFLAGS', 'CFLAGS', 'LDFLAGS']:
        if opt_key in env_dict:
            v = env_dict[opt_key]
            o.write('%s="%s"\nexport %s\n' % (opt_key, v, opt_key))
    for k in ['BOOST_ROOT', 'BJAM_PARENT', 'PHYCAS_ROOT', 'OSTYPE', 'NCL_INSTALL_DIR', 'NCL_ALREADY_INSTALLED']:
        v = env_dict[k]
        o.write('%s="%s"\nexport %s\n' % (k, v, k))
    if IS_MAC:
        o.write('DYLD_LIBRARY_PATH="${NCL_INSTALL_DIR}/lib/ncl:${PHYCAS_ROOT}/phycas/Conversions"\nexport DYLD_LIBRARY_PATH\n')
    else:
        o.write('LD_LIBRARY_PATH="${NCL_INSTALL_DIR}/lib/ncl:${PHYCAS_ROOT}/phycas/Conversions"\nexport LD_LIBRARY_PATH\n')
    o.write('PATH="${BJAM_PARENT}:${PATH}"\nexport PATH\n')
    o.write('PYTHONPATH="${PHYCAS_ROOT}:${PYTHONPATH}"\nexport PATH')

def_phycas_ver = ('1', '2', '0')
phycas_dir = 'Phycas-' + '.'.join(def_phycas_ver)
phycas_compressed = phycas_dir + '.tar.gz'
PHYCAS_DOWNLOAD_LIST = (CodeResource(name='Phycas',
                                  protocol='git',
                                  url='git@github.com:mtholder/Phycas.git',
                                  uncompressed_name='Phycas',
                                  version=None),
                        CodeResource(name='Phycas',
                                  protocol='git',
                                  url='git://github.com/mtholder/Phycas.git',
                                  uncompressed_name='Phycas',
                                  version=None),
                        CodeResource(name='Phycas',
                                  protocol='http',
                                  url='http://www.eeb.uconn.edu/projects/phycas/downloads/v' + def_phycas_ver[0] + '.' + def_phycas_ver[1] + '/' + phycas_compressed,
                                  compressed_name=phycas_compressed,
                                  uncompress_fn=untar,
                                  uncompressed_name=phycas_dir,
                                  version=def_phycas_ver),)

PHYCAS_DEPENDENCY= Dependency(name='Phycas', 
                             env_var='PHYCAS_ROOT', 
                             search_dir_prefix='phycas',
                             resource_list=PHYCAS_DOWNLOAD_LIST)


def_ncl_ver = ('2', '1', '14')
ncl_dir = 'ncl-' + '.'.join(def_ncl_ver)
ncl_compressed = ncl_dir + '.tar.gz'
NCL_DOWNLOAD_LIST = (CodeResource(name='ncl',
                                  protocol='svn',
                                  url='https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1',
                                  uncompressed_name='ncl-2.1',
                                  version=('2', '1')),
                     CodeResource(name='ncl',
                                  protocol='http',
                                  url='http://sourceforge.net/projects/ncl/files/NCL/' + ncl_dir + '/' + ncl_compressed + '/download',
                                  compressed_name=ncl_compressed,
                                  uncompress_fn=untar,
                                  uncompressed_name=ncl_dir,
                                  version=def_ncl_ver),)

NCL_DEPENDENCY= Dependency(name='NCL', 
                             env_var='NCL_INSTALL_DIR', 
                             search_dir_prefix='ncl',
                             resource_list=NCL_DOWNLOAD_LIST)

def_boost_ver = ('1', '45', '0')
boost_dir = 'boost_' + '_'.join(def_boost_ver)
boost_compressed = boost_dir + '.tar.bz2'
BOOST_DOWNLOAD_LIST = (CodeResource(name='boost',
                                    protocol='http',
                                    url='http://sourceforge.net/projects/boost/files/boost/' + '.'.join(def_boost_ver) + '/' + boost_compressed + '/download',
                                    compressed_name=boost_compressed,
                                    uncompress_fn=untar,
                                    uncompressed_name=boost_dir,
                                    version=def_boost_ver),)

BOOST_DEPENDENCY= Dependency(name='boost', 
                             env_var='BOOST_ROOT', 
                             search_dir_prefix='boost_1_',
                             resource_list=BOOST_DOWNLOAD_LIST)

BOOST_DEPENDENCY.acquire()
boost_root = BOOST_DEPENDENCY.full_path
NCL_DEPENDENCY.acquire()
ncl_top = NCL_DEPENDENCY.full_path
PHYCAS_DEPENDENCY.acquire()
phycas_root = PHYCAS_DEPENDENCY.full_path




def get_bjam_par(bjam_src_dir, die_if_not_found):
    bjam_src_subdir = [i for i in os.listdir(bjam_src_dir) if os.path.isdir(os.path.join(bjam_src_dir, i))]
    if IS_MAC:
        bjam_par_pref = 'bin.macosx'
    else:
        bjam_par_pref = 'bin.linux'
    poss_bjam_par_list = [d for d in bjam_src_subdir if d.lower().startswith(bjam_par_pref)]
    if len(poss_bjam_par_list) == 0:
        if die_if_not_found:
            die("Expecting to find a directory starting with", bjam_par_pref, "in", bjam_src_dir, " The build of bjam may have failed")
        
    elif len(poss_bjam_par_list) > 1:
        if die_if_not_found:
            die("Multiple directories starting with", bjam_par_pref, "found in", bjam_src_dir, " Remove (or rename the ones that you do not want to use")
    else:
        bjam_par = poss_bjam_par_list[0]
        if os.path.exists(os.path.join(bjam_src_dir, bjam_par, 'bjam')):
            return os.path.join(bjam_src_dir, bjam_par)
        if die_if_not_found:
            die("bjam not found in ", bjam_par, " The build of bjam may have failed")

    return None

def find_unique_name(fn):
    if os.path.exists(fn):
        n = 1
        while True:
            fp = fn + str(n)
            if not os.path.exists(fp):
                return fp
            n += 1
    return fn

bjam_src_dir = os.path.join(boost_root, 'tools', 'build', 'v2', 'engine', 'src')
bjam_par  = get_bjam_par(bjam_src_dir, False)
if not bjam_par:
    if not do_system(['sh', 'build.sh'], dir=bjam_src_dir):
        die("Could not build bjam")
    bjam_par  = get_bjam_par(bjam_src_dir, True)

ncl_build_dir = os.path.join(ncl_top, 'build')
ncl_install_dir = os.path.join(ncl_build_dir, 'installed')
if not os.path.isdir(os.path.join(ncl_install_dir, 'include', 'ncl', 'ncl.h')):
    if not os.path.exists(os.path.join(ncl_top, 'configure')):
        if not do_system(['sh', 'bootstrap.sh'], dir=ncl_top):
            die("Could not build NCL's configure using bootstrap.sh")
    cfgcommand_fn = os.path.abspath(os.path.join(ncl_build_dir, 'cfgcommand.sh'))
    if not os.path.exists(cfgcommand_fn):
        if not os.path.exists(ncl_build_dir):
            os.mkdir(ncl_build_dir)
        cfgcommand_fo = open(cfgcommand_fn, 'w')
        cfgcommand_fo.write('''#!/bin/sh
../configure --prefix=`pwd`/installed --disable-shared
''')
        cfgcommand_fo.close()
    if not do_system(['sh', 'cfgcommand.sh'], dir=ncl_build_dir):
        die("Could not configure NCL")
    if not do_system(['make', '-j' + NUM_PROC], dir=ncl_build_dir):
        die("Could not build NCL")
    if not do_system(['make', 'install'], dir=ncl_build_dir):
        die("Could not install NCL")
    
    



env = dict(os.environ)
env['BOOST_ROOT'] = boost_root
env['BJAM_PARENT'] = bjam_par
env['NCL_INSTALL_DIR'] = ncl_install_dir
env['NCL_ALREADY_INSTALLED'] = 'true'
if IS_MAC:
    env['OSTYPE'] = 'darwin'
else:
    env['OSTYPE'] = 'linux'
    
env['PHYCAS_ROOT'] = phycas_root

env_fn = find_unique_name('phycas_env.sh')
status('Writing', env_fn)
o = open(env_fn, 'w')
write_sourceable_file(o, env)
o.close()
status('You should be able to "source" ', env_fn, " and then build Phycas by cd'ing to", phycas_root, " and invoking bjam")

build_fn = find_unique_name('build_phycas.sh')
status('Writing', build_fn)
o = open(build_fn, 'w')
o.write('''#!/bin/sh
set -x 
source '%s'
cd '%s'
bjam -j%s release
''' % (os.path.abspath(env_fn), os.path.abspath(phycas_root), NUM_PROC))
o.close()
status('You should be able to rebuild by phycas by running ', build_fn)

