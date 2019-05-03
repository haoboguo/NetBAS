#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cp ../david.id.csh .
    ./david.id.csh
    cd ..
end
