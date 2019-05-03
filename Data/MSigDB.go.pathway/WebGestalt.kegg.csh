#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../WebGestalt.kegg.R | sed "s/changeme/$list[$i]/g" > WebGestalt.kegg.R
    R -f WebGestalt.kegg.R
    cd ..
end
