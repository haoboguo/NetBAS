#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cp ../Venn.Diagram.R .
    R -f Venn.Diagram.R
    cd ..
end
