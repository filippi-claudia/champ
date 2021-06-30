#! /bin/bash
EXECUTABLE=$1
FILE=$2
echo executable $EXECUTABLE 
echo file $FILE 

$EXECUTABLE < $FILE
