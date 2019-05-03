#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat $list[$i].csv >> ../all.proteins.csv
    cd ..
end
