#!/bin/csh -f

set NetFile = $1
set DavFile = $2
set list = `cat $NetFile | sed 's/"//g' | awk '{print $1}'`
set i = 1
while ($i < 14)
  @ i++
  set pval = `grep "$list[$i]" $DavFile | sed 's/ /_/g' | awk '{print $5}'`
  echo $list[$i] $pval
end
