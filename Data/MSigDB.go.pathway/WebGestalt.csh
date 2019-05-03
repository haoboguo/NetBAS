#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../WebGestalt.R | sed "s/changeme/$list[$i]/g" > WebGestalt.R
    R -f WebGestalt.R
    cd ..
end
