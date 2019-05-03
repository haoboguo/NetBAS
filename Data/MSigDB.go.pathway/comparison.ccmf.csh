#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cd $list[$i]
    cp ../gocheck.csh .
    ./gocheck.csh HMK.cc.new.txt DAVID.cc.all.txt
    cp ../comparison.cc.R .
    R -f comparison.cc.R

    ./gocheck.csh HMK.mf.new.txt DAVID.mf.all.txt
    cp ../comparison.mf.R .
    R -f comparison.mf.R

    cd ..
end
