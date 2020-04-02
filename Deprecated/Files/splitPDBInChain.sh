#!/bin/bash

for file in ./PDBs/*
do
  filename=$(basename -- "$file")
  extension="${filename##*.}"
  filename="${filename%.*}"

  cp $file ./PDBs/Chains/${filename}
  mkdir -p "./PDBs/Chains/${filename}"
  perl PDBtoSplitChain.pl -i "$file" -o "./PDBs/Chains/${filename}/${filename}"
done