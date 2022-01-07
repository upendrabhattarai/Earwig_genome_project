# 2. Concatenating, filtering, and classifying repeats

Merging individual repeat library from `RepeatModuler`, `LTRharvest&LTRdigest`, `TransposonPSI`, `SINE database`

```
cat path/to/repeatmoduler/output/consensi.fa.classified \
    path/to/LTRharvest&LTRdigest/pipeline/output/EW_assembly_ltrh.sorted.ltrd.filtered.sequences.fasta \
    path/to/TransposonPSI/output/EW_assembly_renamed.TPSI.allHits.chains.bestPerLocus.fa \
    path/to/SINE/database/SINEs.bnk > Combined.rep.library.fasta
```
