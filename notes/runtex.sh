#!/bin/bash

if [ $# -eq 0 ]
then
    echo usage: 
    echo "      "runtex.sh jobname"          "runs % pdflatex -output-directory files jobname.tex
    echo "      "runtex.sh jobname bib"      "runs % pdflatex -output-directory files jobname.tex
    echo "                                   "then % bibtex jobname
    echo "                                   "followed by % pdflatex -output-directory files jobname.tex twice
elif [ "$2" = "bib" ]
then
    pdflatex -output-directory files $1.tex
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    if [ -e "files/$1Notes.bib" ] 
    then
	cp files/$1Notes.bib . 
    fi
    bibtex files/$1
    if [ -e "$1Notes.bib" ] 
    then
	rm $1Notes.bib
    fi
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    pdflatex -output-directory files $1.tex
    echo
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo
    pdflatex -output-directory files $1.tex
else
    pdflatex -output-directory files $1.tex
fi

if [ -e "files/$1.pdf" ] 
then
    mv files/$1.pdf .
fi

