#!/bin/bash
find . -xtype l 2> tmp0.txt
cat tmp0.txt | awk '{print $2}' | sed 's/`*//g' | sed 's/..$//' >> tmp1.txt
for i in $(cat tmp1.txt); do unlink $i; done
rm tmp0.txt tmp1.txt
