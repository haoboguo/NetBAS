#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../circular.R | sed "s/changeme/$list[$i]/g" > circular.R
    R -f circular.R
    cd ..
end
