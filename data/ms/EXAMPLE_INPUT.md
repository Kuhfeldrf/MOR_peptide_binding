# Expected mass-spec input schema

The Stage 1 parser reads a TSV with at least the columns below. Extra columns
are carried through into the metadata table untouched.

| Column        | Type   | Required | Notes                                    |
|---------------|--------|----------|------------------------------------------|
| `peptide_id`  | string | yes      | Unique within the file.                  |
| `sequence`    | string | yes      | One-letter code, canonical residues only.|
| `id_score`    | float  | no       | Instrument/search-engine confidence.     |
| `charge`      | int    | no       | Precursor charge state.                  |
| `rt_min`      | float  | no       | Retention time in minutes.               |
| `source`      | string | no       | Sample or fraction identifier.           |

## Status

**Real peptide data is in use.** This schema describes the generic MS input
contract, but Stage 1 does not depend on a synthetic stand-in: the pilot reads
the curated peptide set carried over from the prior ORCA work. See
`data/reference/README.md` for the source and its provenance.

`example_peptides.tsv` is retained only as a three-line format illustration for
anyone wiring up a new instrument export. It is labelled `SYNTHETIC_EXAMPLE` in
its `source` column and is not read by any rule.
