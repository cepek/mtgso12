#!/bin/sh

PROJECT_NAME=randm-a

printf 'Generating randmat test data and plot. Do you want to continue? (y/n) '
read answer

if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo Yes ... it will take some time
    echo
else
    echo No
    exit 0;
fi


if [ ! -f ../../build/randmat ];
then
    echo "../../build/randmat could not be found"
    exit
fi


echo Cleaning data files ...
rm -vf *$PROJECT_NAME.txt *$PROJECT_NAME.dat $PROJECT_NAME-files.txt

echo
echo Running octave ...
octave < $PROJECT_NAME.m
echo

echo
echo Running tests ...

for i in *$PROJECT_NAME.dat
do
    sed '/^#/d' $i > ${i%%.dat}.txt
done

ls *$PROJECT_NAME.txt > $PROJECT_NAME-files.txt
../../build/randmat . $PROJECT_NAME-files.txt

gnuplot plot.gp
#evince $PROJECT_NAME.eps &
