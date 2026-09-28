# Contributing

This repository is built from Sapote-Mamey by `build_rggmci_package.py`. The scorer (`core.py`) and the
lifted helpers (`_parsers.py`, `_ids.py`, `_blastp_io.py`) are copies of Sapote-Mamey source. Don't edit them
here: a change made here would be lost on the next build and would make this package disagree with the
engine. Report problems in them as issues; fixes go into Sapote-Mamey and reach this repository on the next
build.

The reader, grouping, BLASTp round trip, command line, tests and docs are package-only and can be changed
here, with a matching change to the templates in Sapote-Mamey.

`tests/test_example_matches_engine.py` pins the package to the engine's own output on the public example, so
any drift fails the tests.
