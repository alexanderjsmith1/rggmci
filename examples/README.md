# Example

`VWPH00000000.1_subset.antismash.zip` holds five regions, with their ClusterBlast files, cut from the public antiSMASH result for NCBI WGS VWPH00000000.1 (*Saccharopolyspora*). It gives 1 scored pairs, 0 of them HIGH.

```bash
rggmci examples/VWPH00000000.1_subset.antismash.zip --out result.json --pairs pairs.tsv
```

`VWPH00000000.1_subset_expected_engine.json` is the Sapote-Mamey engine's own result on this file; `tests/test_example_matches_engine.py` checks that the package reproduces it exactly.
