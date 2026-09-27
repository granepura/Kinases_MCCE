# Experimental pKa values used in Table S6

Provenance for the measured pKa values compared against the calculated solution
protonation states in Table S6 of *Kinase inhibitors can change protonation or
tautomeric state upon binding*.

## Why this directory exists

The calculations in this work run at a single pH and determine no pKa. What Table S6
compares is the **solution reference** each inhibitor is assigned, which is a Boltzmann
distribution over Epik state penalties, against pKa values measured in the literature.
Every measured value therefore comes from outside this project, and this directory records
exactly where each one came from so the comparison can be checked without repeating the
search.

## Files

- `pka_provenance.csv` — one row per value, including the negative results.
- `source_documents/` — the primary documents, named `<inhibitor>_<BRAND>_<source>_<year>.pdf`.
- `verify_pka_provenance.py` — re-checks every row against its PDF.

## Manifest columns

| Column | Meaning |
|---|---|
| `inhibitor` | Inhibitor name as used in the manuscript. |
| `site_as_reported` | The ionizable group, and whether the **source** assigned it or we assigned it by inspection of the structure. Six sources give a bare number with no site. |
| `pka_calc_used` | The calculated value in Table S6, for cross-reference. |
| `pka_exp_used` | The measured value used in Table S6. |
| `pka_as_reported` | The value as the source states it, including salt form, temperature and ionic strength where given. |
| `page` | Page of the PDF carrying the sentence. |
| `quoted_sentence` | The sentence itself, verbatim. For rows with no value this states the negative result. |
| `source_document` | Document title. |
| `url` | Where it was retrieved. |
| `local_pdf` | Filename in `source_documents/`. |
| `sha256` | Checksum of that file, so a re-download can be shown to be the same bytes. |
| `verification` | Result of the automated check. |
| `accessed` | Retrieval date. |

## Verifying

```bash
python3 verify_pka_provenance.py
```

Recomputes every checksum, confirms each quoted sentence is present in its PDF on the
stated page, and confirms that the documents recorded as carrying no pKa really carry
none. Requires PyMuPDF.

## What is and is not here

Twelve of the eighteen inhibitors have a measured value. Six do not: bosutinib,
cobimetinib, lapatinib, osimertinib, sorafenib and tofacitinib. For four of those the
absence is documented here as a negative result rather than merely asserted, because the
paper states that the coverage is partial. Cobimetinib and osimertinib were searched
without success and have no document recorded.

Two caveats travel with these numbers. The lenvatinib and sunitinib values are reported
for the mesylate and malate salts, and the nilotinib pKa2 is reported by its source as
estimated rather than determined. Both are noted in the Table S6 footnote.

The imatinib value is the only one from the primary literature rather than a regulatory
document, and it is the **one row with no local copy and no checksum**: the publisher
returned HTTP 403, so the value was taken from the abstract record and still needs
checking against the article PDF.

## A note on the PDFs

These are third-party regulatory and manufacturer documents, reproduced here unmodified
so the quoted sentences can be checked. They are not our work and carry their original
terms. The manifest alone, with its URLs and checksums, is sufficient to verify every
value; if redistribution is ever a concern, deleting `source_documents/` leaves the
provenance record intact and `verify_pka_provenance.py` will report the files as missing
rather than failing.
