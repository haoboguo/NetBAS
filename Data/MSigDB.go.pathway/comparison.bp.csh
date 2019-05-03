#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cp ../gocheck.csh .
    ./gocheck.csh HMK.bp.new.txt DAVID.bp.all.txt
    cp ../comparison.bp.R .
    R -f comparison.bp.R
    cd ..
end
