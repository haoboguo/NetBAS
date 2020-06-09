#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    set num =  `cat $list[$i]/$list[$i].csv | wc -l`
    set numb = `echo $num - 1 | bc -l`
    echo $list[$i] $numb >> gene.count.csv
end
