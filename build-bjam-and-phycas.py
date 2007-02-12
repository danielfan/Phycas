#! /usr/bin/env python
import os
import sys
import subprocess

bundled_boost = "boost_1_33_1"
phycas_version_str = "0.11.0"
bundled_phycas = "Phycas-" + phycas_version_str

################################################################################
# Some globals that control script behavior
################################################################################
interactive_mode = True
n_prompt_tries = 10
default_verbosity = False
is_windows = sys.platform.upper().startswith('WIN')
is_mac = sys.platform.upper() == 'DARWIN'
_parent_dir, _prog = os.path.split(sys.argv[0])
_parent_dir = os.path.abspath(_parent_dir)

################################################################################
# Minimal logging functions
# NOTE: calling error reports an error and exits the script!
################################################################################
def info_message(s):
    "Prints message with prefix for clarifying the source of the message."
    print _prog + ': ' + s

def error(s):
    """Displays the message (with program name and ERROR tag) then exits with 
    error code."""
    sys.exit(_prog + ': ERROR\n' + s + '\n' + _prog + ': Exiting due to ERROR.')

def warn(s):
    "Displays the message (with program name and WARNING tag)"
    sys.stdout.write(_prog + ': WARNING\n' + s + '\n')

def printEnv(environ):
    "Prints keys and values in environ in an easy-to-read format"
    for k, v in environ.iteritems():
        info_message('%-20s = %s' % (k,v))


def versionTupleCmp(f, s):
    """__cmp__ function for version tuples.
    cannot use lexigraphic cmp because:
        1.1 < 1.04
    can't convert the whole thing to a decimal because
        1.10 > 1.1.11"""
    if len(f) == 0:
        return len(s) != 0 and -1 or 0
    for ind, v in enumerate(f):
        if ind+1 > len(s):
            return 1
        r = cmp(v, s[ind])
        if r != 0:
            return r
    if len(s) > len(f):
        return -1
    return 0

def set_env(k, val, env_dict, sec_dict=None):
    if k in env_dict:
        old = env_dict[k]
        info_message("Resetting %s from %s to %s" % (k, old, val))
    else:
        info_message("Setting %s %s" % (k, val))
    env_dict[k] = val
    if sec_dict is not None:
        sec_dict[k] = val

def verify_path_exists(p):
    if not os.path.exists(p):
        error('Required path "%s" not found.' % p)
        return False
    return True

def check_for_bjam(boost_root, required=False):
    former_dir = os.curdir
    bjam_src_par = os.path.join(boost_root, "tools", "build", "jam_src")
    verify_path_exists(bjam_src_par)
    os.chdir(bjam_src_par)
    try:
        sub_paths = os.listdir(os.curdir)
        for p in sub_paths:
            if p.startswith("bin.") and os.path.isdir(p):
                s = os.path.join(bjam_src_par, p, "bjam")
                if os.path.exists(s):
                    a = ["$BOOST_ROOT", "tools", "build", "jam_src", p, "bjam"]
                    return os.path.join(*a)
        if required:
            error("bjam not found in a subdirectory of %s" % bjam_src_par)
        return None
    finally:
        os.chdir(former_dir)

def build_bjam(boost_root):
    former_dir = os.curdir
    bjam_src_par = os.path.join(boost_root, "tools", "build", "jam_src")
    verify_path_exists(bjam_src_par)
    os.chdir(bjam_src_par)
    try:
        if subprocess.call(["/bin/sh", "build.sh"]) != 0:
            error("Could not build bjam from %s" % bjam_src_par)
        return check_for_bjam(boost_root, True)
    finally:
        os.chdir(former_dir)

def write_phycas_build_sh(phypy_dir, bjam_path, env):
    former_dir = os.curdir
    expanded_phypy_dir = os.path.expandvars(phypy_dir)
    os.chdir(expanded_phypy_dir)
    try:
        f = file("build_env.sh", "w")
        f.write("#!/bin/sh\n\n")
        ks = ["PYTHON_ROOT", "PYTHON_VERSION", "BOOST_ROOT", "BUILD", "TOOLS", "PHYCAS_ROOT"]
        for k in ks:
            v = env[k]
            f.write('%s="%s"\nexport %s\n' % (k, v, k))
        bjam_par = os.path.split(bjam_path)[0]
        f.write('\nPATH="%s%s${PATH}"\nexport PATH\n\n' % (bjam_par, os.pathsep))
        f.close()

        f = file("build.sh", "w")
        f.write("#!/bin/sh\n\n")
        k = "PHYCAS_ROOT"
        v = env[k]
        f.write('%s="%s"\nexport %s\n' % (k, v, k))
        f.write('\nsource "%s" || exit\n' % os.path.join(phypy_dir, "build_env.sh"))
        f.write('\nexec python "%s"\n' % os.path.join(phypy_dir, "dojam.py"))

        subprocess.call(["chmod", "775", "build.sh"])
        return os.path.join(phypy_dir, "build.sh")
    finally:
        os.chdir(former_dir)
    
if __name__ == '__main__':
    ld_lib_path_var = is_mac and "DYLD_LIBRARY_PATH" or "LD_LIBRARY_PATH"

    starting_dir = os.curdir
    os.chdir(_parent_dir)

    env = os.environ
    init_pythonpath = env.get("PYTHONPATH")
    init_ld_lib_path = env.get(ld_lib_path_var)
    added_to_env = {}
    boost_root = env.get("BOOST_ROOT")
    if boost_root:
        info_message("BOOST_ROOT set to %s" % boost_root)
    else:
        boost_root = os.path.abspath(os.path.join(_parent_dir, bundled_boost))
        set_env("BOOST_ROOT", boost_root, env, added_to_env)

    bjam_path = check_for_bjam(boost_root)
    if not bjam_path:
        bjam_path = build_bjam(boost_root)
    expanded_bjam = os.path.expandvars(bjam_path)
    set_env("PYTHON_ROOT", sys.prefix, env, added_to_env)
    v = sys.version_info[:2]
    python_version = str(v[0]) + '.' + str(v[1])
    set_env("PYTHON_VERSION", python_version, env, added_to_env)
    
    
    build_mode = "release"
    set_env("BUILD", build_mode, env, added_to_env)
    tools = is_mac and "darwin" or "gcc"
    set_env("TOOLS", tools, env, added_to_env)
    
    phycas_root = env.get("PHYCAS_ROOT")
    if phycas_root:
        info_message("PHYCAS_ROOT set to %s" % phycas_root)
    else:
        phycas_root = os.path.join(_parent_dir, bundled_phycas)
        set_env("PHYCAS_ROOT", phycas_root, env, added_to_env)
    phypy_dir = os.path.join("$PHYCAS_ROOT", "phypy")
    expanded_phypy_dir = os.path.expandvars(phypy_dir)
    
    
    # write the build script
    phypy_build_sh = os.path.join(expanded_phypy_dir, "build.sh")
    if not os.path.exists(phypy_build_sh):
        write_phycas_build_sh(phypy_dir, bjam_path, env)
    
    # do the build
    if subprocess.call(["/bin/sh", phypy_build_sh]) != 0:
        error("Could not build using %s" % phypy_build_sh)
    
    
    
    #write the env scripts
    os.chdir(_parent_dir)
    bash_vars = True
    shell_var = os.path.split(env["SHELL"])[-1]
    if shell_var == "csh" or shell_var == "tcsh":
        bash_vars = False

    # find the libboost_python library that will have to be loaded 
    #  dynamically
    lb_ext = is_mac and "dylib" or "so"
    lb_name = "libboost_python." + lb_ext
    p = [phypy_dir, 
         "bin", 
         "boost", 
         "libs", 
         "python", 
         "build", 
         lb_name, 
         tools, 
         build_mode,
         "shared-linkable-true",
         ]
    libboost_par = os.path.join(*p)
    exp_libboost_par = os.path.expandvars(libboost_par)
    libboost_path = os.path.join(exp_libboost_par, lb_name)
    if not os.path.exists(libboost_path):
        m = "%s not found load library path not determined!" % libboost_path
        sys.exit(m)

    suffix = bash_vars and ".sh" or ".csh"
    to_source = "phycas_env" + suffix
    test_file = "test_phycas" + suffix
    rc_file = bash_vars and ".profile" or ".cshrc"
    source_cmd = 'source "%s"' % os.path.abspath(to_source)
    tests_dir = os.path.join(expanded_phypy_dir, "phypy", "Tests")
    test_cmd = """cd %s
    python runall.py
    python doctestall.py
    """ % tests_dir
    print """
To use Phycas, you will have to add the update the PYTHONPATH and variables 
in your environment. You can do this by "sourcing" the file %s
using the command:

    %s

Note: you may also want to add this source command to your %s 
file so that
the file is read and the environment is configured whenever you log in.
""" % ( to_source, 
        source_cmd,
        os.path.expandvars(os.path.join("$HOME", rc_file)),
       )

    print """
IMPORTANT:
    If you move the Phycas directory, you will have to update the %s file
    to reflect the new location of Phycas"
    
After sourcing %s you should be able to run the following command without
errors:
    python -c "import phypy"

The more rigorous, but time-consuming tests can be invoked with:
    
    %s

""" % ( to_source,
        tests_dir,
        test_cmd,
      )
    
    old_ld_lib_path = "%s${%s}" % (os.pathsep, ld_lib_path_var)
    full_ld_lib_path = "%s%s" % (exp_libboost_par, old_ld_lib_path)
    if bash_vars:
        lib_setting = '%s="%s"; export %s' % (ld_lib_path_var,
                                    full_ld_lib_path,
                                    ld_lib_path_var)
    else:
        lib_setting = 'setenv %s "%s"' % (ld_lib_path_var, full_ld_lib_path)
    print lib_setting

    otherpp = os.pathsep + "${PYTHONPATH}"
    fullpp = "%s%s" % (expanded_phypy_dir, otherpp)
    if bash_vars:
        pypath_setting = 'PYTHONPATH="%s"; export PYTHONPATH' % fullpp
    else:
        pypath_setting = 'setenv PYTHONPATH "%s"' % fullpp
    print pypath_setting
    
    if os.path.exists(to_source):
        info_message("%s found.  It will not be replaced." % to_source)
    else:
        f = open(to_source, "w")
        f.write("#! /bin/sh\n%s\n%s\n" % (pypath_setting, lib_setting))
        f.close()
        info_message("%s written." % to_source)

    if os.path.exists(test_file):
        info_message("%s found.  It will not be replaced." % test_file)
    else:
        f = open(test_file, "w")
        f.write("#! /bin/sh\n%s\n" % test_cmd)
        f.close()
        info_message("%s written." % test_file)
