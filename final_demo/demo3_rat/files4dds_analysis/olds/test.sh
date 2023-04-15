#!/bin/bash

line='rn6.refFlat.txt'
replace='W' 
finds='6'
echo $line | sed -e "s/6/${replace}/g"

result=${line//$finds/$replace}
echo ${result}
