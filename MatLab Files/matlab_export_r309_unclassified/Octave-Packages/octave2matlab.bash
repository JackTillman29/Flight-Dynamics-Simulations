#!/bin/bash

sed --in-place *.m -e 's,##,%,g' \
                   -e 's,#,%,g' \
                   -e 's,!,~,g' \
                   -e 's,endif,end,g' \
                   -e 's,endfor,end,g' \
                   -e 's,endswitch,end,g' \
                   -e 's,endfunction,end,g' \
                   -e 's,__,,g'



for i in `ls __*__.m`;
do
  k=${i%__*}   #remove suffix underscores
  k=${k#*__}   #remove prefix underscores
  k=$k.m       #add .m file extension back in
  echo "renaming $i to $k"
  mv $i $k
done

for i in `ls __*.m`;
do
  k=${i#*__}   #remove prefix underscores
  echo "renaming $i to $k"
  mv $i $k
done
