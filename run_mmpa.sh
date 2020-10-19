#!/usr/bin/bash
set -e 

jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --inplace --execute extract_scaffold_memory.part1.ipynb
jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --inplace --execute extract_scaffold_memory.part2.ipynb

python mmpa_fragmentation.py
parallel < to_fragment/jobs 

python mmpa_indexing.py
parallel < fragmented_smiles/jobs

jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --inplace --execute analyze_mmpa.part1.ipynb
jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --inplace --execute analyze_mmpa.part2.ipynb

