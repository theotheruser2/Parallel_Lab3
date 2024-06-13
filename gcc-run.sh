#!/bin/bash

touch ./res-gcc/parallel3.txt
touch ./res-gcc/seq3.txt
for (( c=10; c<=400; c+=10 ))
do  
   ./lab3 $c >> ./res-gcc/parallel3.txt
   ./lab3seq $c >> ./res-gcc/seq3.txt

done

