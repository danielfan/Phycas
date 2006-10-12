#!python

import sys,os,distutils.file_util

if len(sys.argv) > 0 and sys.argv[1] == '-install':
    site_packages = os.path.join(sys.prefix , 'Lib', 'site-packages')

    # Create a directory to hold Phycas-related shortcuts
    common_programs = get_special_folder_path('CSIDL_COMMON_PROGRAMS')
    Phycas_shortcuts = os.path.join(common_programs, 'Phycas')
    if not os.path.isdir(Phycas_shortcuts):
        os.mkdir(Phycas_shortcuts)
        directory_created(Phycas_shortcuts)

    # Create shortcut to runall.bat file
    target = os.path.join(site_packages, 'phypy', 'Tests', 'runall.bat')
    description = 'Batch file that tests Phycas'
    filename = 'Test Phycas.lnk'
    arguments = ''
    workdir = os.path.join(site_packages, 'phypy', 'Tests')
    create_shortcut(target, description, filename, arguments, workdir)
    dest_name, ok = distutils.file_util.copy_file(filename, Phycas_shortcuts)
    file_created(os.path.join(Phycas_shortcuts, filename))
    os.remove(filename)

    # Create shortcut to uninstall
    target = '"' + os.path.join(sys.prefix, 'RemovePhycas.exe') + '"'
    description = 'This shortcut will uninstall Phycas from your computer'
    filename = 'Uninstall Phycas.lnk'
    arguments = '-u "' + os.path.join(sys.prefix, 'Phycas-wininst.log') + '"'
    workdir = sys.prefix
    create_shortcut(target, description, filename, arguments, workdir)
    dest_name, ok = distutils.file_util.copy_file(filename, Phycas_shortcuts)
    file_created(os.path.join(Phycas_shortcuts, filename))
    os.remove(filename)

    print 'Please use the "Test Phycas" shortcut and let us know of any problems reported'
elif len(sys.argv) > 0 and sys.argv[1] == '-remove':
    pass
else:
    print "expecting either '-install' or '-remove' to be supplied as an argument"
