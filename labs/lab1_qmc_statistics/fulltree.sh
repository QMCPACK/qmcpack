#!/bin/sh
#############################################
# Script that displays a recursive formatted folder and file listing
# @author Corbin
# @site iamcorbin.net
#Folder Seperator
BREAK='-------------------------------------------------------------------------------------'

#Optional: if a folder is passed as an argument, run fulltree on that folder rather than the current folder
if [ "$1" != "" ]
   then cd "$1"
   fi
pwd

## Recursive Directory Listing with files
 # 1- preserve directories from being removed in 2 & 3
 # 2- strip first 4 columns
 # 3- strip size and date
 # 4- prepend '  -- ' on each line
 # 5- remove '  -- ' from directories
 # 6- remove extra lines
 # 7- Insert a line break after directories
 # 8- Put a | at the beginning of all lines
 # 9- Indent and process 1st level sub dirs
 #10- Indent and process 2nd level sub dirs
ls -Rhl | sed \
    -e 's/^\.\//x x x x 00:00 |-/' \
    -e 's/^\([^\ ]*.\)\{4\}//' \
    -e 's/.*[0-9]\{2\}:[0-9]\{2\}//' \
    -e 's/^/  -- /' \
    -e 's/\ \ --\ \ |-//'  \
    -e '/--\ $/ d' \
    -e '/^[^ ]/ i\'$BREAK \
    -e 's/^/| /' \
| sed -e '/[^/]*\//,/'$BREAK'/ s/^|/\t&/' -e '/^\t/,/'$BREAK'/ s/'$BREAK'/\t&/' -e 's/[^/]*\//\t\| /' \
| sed -e '/[^/]*\//,/'$BREAK'/ s/^\t|/\t&/' -e '/^\t\t/,/'$BREAK'/  s/'$BREAK'/\t&/' -e 's/[^/]*\//\t\t\| /' \
| sed -e '/[^/]*\//,/'$BREAK'/ s/^\t\t/\t&/' -e 's/[^/]*\//\t\t\t\| /'
echo $BREAK
