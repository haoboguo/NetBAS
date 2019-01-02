#!/bin/csh -f

set list = `cat list`
set i = 0
while ($i < $#list)
    @ i++
    cat $list[$i]/HMK.kegg.enriched.txt | sed "s/Z-score/zscore/g" > $list[$i]/HMK.kegg.full.csv
    cat $list[$i]/HMK.bp.full.txt | sed "s/Z-score/zscore/g" > $list[$i]/HMK.bp.full.csv
    cat $list[$i]/HMK.cc.full.txt | sed "s/Z-score/zscore/g" > $list[$i]/HMK.cc.full.csv
    cat $list[$i]/HMK.mf.full.txt | sed "s/Z-score/zscore/g" > $list[$i]/HMK.mf.full.csv
end
