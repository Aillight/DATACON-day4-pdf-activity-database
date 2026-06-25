# Validation Report

## Summary

- Source PDFs: 8
- Final measurement records: 99
- Manual review records: 59
- Duplicate candidate groups: 5
- Conflict groups: 2
- Records without resolved structure identifiers: 99

## Top Validation Issues

- multiple_possible_targets: 32
- compound_not_resolved: 20
- duplicate_candidate: 10
- unit_not_resolved: 10
- multi_value_respectively_context: 9
- reference_compound_context: 6
- potential_conflict_group: 6
- target_not_resolved: 2

## Records By Source

- 1-s2.0-S0960894X09006258-main.pdf: 37
- 1-s2.0-S0960894X09003679-main.pdf: 28
- 1-s2.0-S0223523412000402-main.pdf: 14
- 1-s2.0-S0960894X09006222-main.pdf: 11
- design-synthesis-and-biological-evaluation-of-6α-and-6β-n-heterocyclic-substituted-naltrexamine-derivatives-as-μ-opioid.pdf: 5
- further-studies-at-neuropeptide-s-position-5-discovery-of-novel-neuropeptide-s-receptor-antagonists.pdf: 2
- 1-s2.0-S0960894X08008214-main.pdf: 1
- synthesis-and-evaluation-of-dibenzothiazepines-a-novel-class-of-selective-cannabinoid-1-receptor-inverse-agonists.pdf: 1

## Records By Target

- NOP/ORL1: 39
- MOR: 12
- hERG: 9
- MOR;DOR: 7
- NOP/ORL1;MOR: 5
- NOP/ORL1;DOR;KOR: 5
- DOR;KOR: 5
- KOR;MOR: 4
- KOR: 4
- MOR;DOR;KOR: 3
- unresolved: 2
- DOR;KOR;NOP/ORL1: 1
- MOR;NOP/ORL1: 1
- KOR;MOR;DOR: 1
- CB1: 1

## Extraction Methods

- pypdf_text_regex_v1: 62
- pypdf_orl1_table_regex_v1: 37

## Conflict Groups

- conflict_6805676f9f: 16 / NOP/ORL1 / ED50 -> =0.3 mg/kg, =3.2 mg/kg, =8.8 mg/kg
- conflict_9012ee57fc: 15 / NOP/ORL1 / ED50 -> =0.02 mg/kg, =0.3 mg/kg, =6.7 mg/kg

## Notes

- Values are not averaged. Potential duplicates and conflicts remain visible for curation.
- `manual_review.csv` is the queue for curator corrections.
- `compound_identity_map.csv` is the offline place to add SMILES, InChI, PubChem CID and CAS RN.
