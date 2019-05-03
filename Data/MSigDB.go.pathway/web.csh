#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cp Project*/*enrichedResult.txt WebGestalt.2.txt
    cd ..
end
