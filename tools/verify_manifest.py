from pathlib import Path
import hashlib,csv
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((root/'CHECKSUMS_SHA256.csv').open(encoding='utf-8')))
bad=[r['file'] for r in rows if not (root/r['file']).is_file() or hashlib.sha256((root/r['file']).read_bytes()).hexdigest()!=r['sha256']]
assert not bad,bad
print('Verified',len(rows),'bundled files')
