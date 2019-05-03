#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cat david.bp.new.txt | sed "s/\~/ /g" | awk '{print $1}' > david.bp.id.txt
    cat david.cc.new.txt | sed "s/\~/ /g" | awk '{print $1}' > david.cc.id.txt
    cat david.mf.new.txt | sed "s/\~/ /g" | awk '{print $1}' > david.mf.id.txt
    cd ..
end
