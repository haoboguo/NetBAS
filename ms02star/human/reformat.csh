#!/bin/csh -f

set i = 0

while ($i < 9)
  @ i++
  cat ms02.$i.csv | sed "s/,/ /g" | awk '{print $1,$2}' | sed 's/ /,/' > /Users/yhw543/Documents/GitHub/NetBAS_Codes-Resourses/data/human.ms02/ms02.$i.csv
end
