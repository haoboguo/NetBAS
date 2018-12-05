#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    R -f comparison.bp.R
    R -f comparison.cc.R
    R -f comparison.mf.R
    cd ..
end
