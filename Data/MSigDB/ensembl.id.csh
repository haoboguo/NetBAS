#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat ../ensembl.id.R | sed "s/changeme/$list[$i]/g" > ensembl.id.R
    R -f ensembl.id.R
    cd ../
end
