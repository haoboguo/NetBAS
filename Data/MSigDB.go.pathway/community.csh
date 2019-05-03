#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../community.R | sed "s/changeme/$list[$i]/g" > community.R
    R -f community.R
    cd ..
end
