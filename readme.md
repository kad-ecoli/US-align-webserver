# US-align-webserver #
Source code for the webserver of US-align

## Installation ##
1. The CGI scripts are known to be compatible with python2,7 up to python3.12.
   You may want to modify the shebang line (first line) of 
   `bin/receive*.cgi`` to point to the correct version of python interpreter

2. run ``output/readme.sh`` to setup the correct permission and context

3. you may want to add ``remove_old_job.sh`` to hourly crontab to remove
   old job files automatically.

## Update ##
US-align executable is compiled from the single-file source code at 
``bin/module/USalign.cpp``. After updating this file, recompile it by
```bash
cd bin/module/
make USalign
```
