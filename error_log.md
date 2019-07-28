# Errors to be fixed

1. A circular import error in `make_all_nmer_substitutions.py` not allowing packages from this to be imported into other scripts. The error lies, likely, with the `codon_table()` module within `make_all_nmer_substitutions.py`. 
