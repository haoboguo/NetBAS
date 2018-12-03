#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cp ../MSigDB/$list[$i]/$list[$i].csv $list[$i]/.
end
