#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../pathfindR.R | sed "s/changeme/$list[$i]/g" > pathfindR.R
#    R -f pathfindR.R
    cd ..
end
