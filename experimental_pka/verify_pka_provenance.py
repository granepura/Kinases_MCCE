#!/usr/bin/env python3
"""Re-verify pka_provenance.csv against the PDFs in source_documents/.

Checks, per row: the checksum still matches, the quoted sentence is present in the
named PDF on the stated page, and rows recorded as carrying no pKa really carry none.
Exits non-zero if anything fails. Requires PyMuPDF.
"""
import csv, hashlib, os, re, sys

try:
    import fitz
except ImportError:
    sys.exit('PyMuPDF is required: pip install pymupdf')

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(HERE, 'source_documents')
norm = lambda s: re.sub(r'\s+', ' ', s).strip()

def real_pka_hits(text):
    """pKa mentions that look like a stated value, not 'PK' pharmacokinetics."""
    return [m for m in re.finditer(r'pK\s*a', text)
            if re.search(r'pKa\s*(?:of|value|1|2|=|\d)', text[m.start():m.start() + 24], re.I)]

def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, 'pka_provenance.csv'))))
    bad = skipped = 0
    for r in rows:
        tag = f"{r['inhibitor']:12} {r['source_document'][:42]:44}"
        pdf = r['local_pdf']
        if not pdf:
            print(f'  skip {tag} no local copy'); skipped += 1; continue
        path = os.path.join(DOCS, pdf)
        if not os.path.exists(path):
            print(f'  skip {tag} file absent'); skipped += 1; continue
        if hashlib.sha256(open(path, 'rb').read()).hexdigest() != r['sha256']:
            print(f'  FAIL {tag} checksum changed'); bad += 1; continue
        doc = fitz.open(path)
        full = norm(' '.join(p.get_text() for p in doc))
        quote = r['quoted_sentence']
        if quote.startswith('NEGATIVE RESULT'):
            hits = real_pka_hits(full)
            if hits:
                print(f'  FAIL {tag} expected no pKa, found {len(hits)}'); bad += 1
            else:
                print(f'  ok   {tag} absence confirmed')
        elif norm(quote) not in full:
            print(f'  FAIL {tag} quoted sentence not found'); bad += 1
        else:
            page = r['page']
            where = ''
            if page:
                i = int(page) - 1
                if not (0 <= i < len(doc) and norm(quote) in norm(doc[i].get_text())):
                    where = f' (WARNING: not on stated page {page})'
            print(f'  ok   {tag} quote found{where}')
    print(f'\n{len(rows)} rows: {len(rows) - bad - skipped} verified, {skipped} skipped, {bad} failed')
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main())
