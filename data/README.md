# Data

The trial data are not included in this repository. They are publicly posted by the trial
investigators:

- Doyle J, Alsan M, Skelley N, Lu Y, Cawley J. Replication Data for: Effect of an Intensive
  Food-as-Medicine Program on Health and Health Care Use: A Randomized Clinical Trial. Harvard
  Dataverse, version 2. doi:10.7910/DVN/LZUTU0
- https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LZUTU0

Download the files and place them here:

```
data/raw/dataverse_files/
  analysis_sample.csv        465 patients x 157 variables (deidentified)
  README.pdf                 variable dictionary (from the investigators)
  01_published_tables.do     the investigators' Stata code (used by verification/stata/)
```

Notes used by the code: dates are Stata day numbers (days since 1960-01-01); one baseline
triglyceride value recorded as 9999999 is treated as missing; weight in kg is recomputed from the
pound columns divided by 2.20462 (the posted kg column used 2.2222).
