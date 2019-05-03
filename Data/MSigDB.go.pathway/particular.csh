#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../particular.go.term.R | sed "s/changeme/$list[$i]/g" > particular.go.term.R
    R -f particular.go.term.R
    cd ..
end
