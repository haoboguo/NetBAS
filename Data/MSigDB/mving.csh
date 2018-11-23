#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    mv $list[$i].csv $list[$i]/.
end
