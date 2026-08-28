#!/bin/bash

# This folder is the temporary folder for writing pdb files by cgi script.
# Permission and context of the folder and cgi scripts must be correctly set:

# First, set the permission:

    chmod a+x bin/receive*.cgi
    chmod a+x bin/module/USalign
    chmod a+x bin/module/cif2pdb
    chmod a+x bin/module/pdb2fasta
    touch     log/bywhom.txt
    chmod 777 log/bywhom.txt
    chmod 777 output/

# Second, set the context, which may not be necessary on some system:

    chcon -t httpd_sys_script_exec_t bin/receive*.cgi
    chcon -t httpd_sys_script_exec_t bin/module/USalign
    chcon -t httpd_sys_script_exec_t bin/module/cif2pdb
    chcon -t httpd_sys_script_exec_t bin/module/pdb2fasta
    chcon -t httpd_sys_rw_content_t  output/
    chcon -t httpd_sys_rw_content_t  log/bywhom.txt

# You can check the context by

    ls -laZ bin/receive*.cgi bin/module/USalign bin/module/cif2pdb bin/module/pdb2fasta
    ls -laZ|grep output

