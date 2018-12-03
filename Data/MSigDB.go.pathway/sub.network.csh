#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../sub.network.R | sed "s/changeme/$list[$i]/g" > sub.network.R
    R -f sub.network.R
    cd ..
end
