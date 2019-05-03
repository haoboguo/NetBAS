#!/bin/csh -f

set NetFile = $1
set DavFile = $2
set list = `cat $NetFile | sed 's/"//g' | awk '{print $1}'`
set i = 1
echo "NetBAS DAVID" > net.dav.txt
while ($i < 11)
  @ i++
  set pval = `grep "$list[$i]" $DavFile | sed 's/ /_/g' | awk '{print $5}'`
  echo $list[$i] $pval >> net.dav.txt
end
