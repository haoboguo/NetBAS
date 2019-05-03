#!/bin/csh -f

set list = (bp mf cc)
set i = 0
while ($i < $#list)
  @ i++
  cat david.$list[$i].new.txt | sed "s/\~/ /g" | awk '{print $1}' > david.$list[$i].id.txt
end
