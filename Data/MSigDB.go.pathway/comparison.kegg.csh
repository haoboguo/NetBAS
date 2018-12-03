#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    R -f comparison.kegg.R
    cd ..
end
