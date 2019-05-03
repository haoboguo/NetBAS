#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../neibor.network.R | sed "s/changeme/$list[$i]/g" > neibor.network.R
    R -f neibor.network.R
    cd ..
end
