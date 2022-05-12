#!/bin/csh

set date = `tail -2 $1`
set runtime = `tail -1 $1`
set days = $runtime[1]
set seconds = $runtime[2]
set curr_date = $date[1],$date[2],$date[3],$date[4],$date[5],$date[6]

cat > modify_file << ENDNL
/current_date/ c\
  current_date = $curr_date
/months / c\
  months   =    0 
/days / c\
  days     =    $days
/hours / c\
  hours    =    0
/minutes / c\
  minutes  =    0
/seconds / c\
  seconds  =    $seconds
ENDNL

sed -f modify_file $2 > file_mod.txt

mv file_mod.txt $2
 
