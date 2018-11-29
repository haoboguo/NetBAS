#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cat ../MSigDB/$list[$i]/$list[$i].geneset.csv | sed 's/"//g' | sed 's/,/ /g' | awk '{print $1}' > $list[$i]/$list[$i].geneset.txt
end
