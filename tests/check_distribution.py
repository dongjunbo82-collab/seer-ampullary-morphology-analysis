"""Fail closed if the code-only distribution contains data or local personal metadata."""
from pathlib import Path
import re,json
root=Path(__file__).resolve().parents[1]
top={'README.md','LICENSE','.gitignore','CITATION.cff','.zenodo.json','run_analysis.R'}
exts={'.R','.md','.py'}
bad=[]
patterns=[r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]',r'(?i)C:[\\/]Users',r'(?i)gh[pousr]_[A-Za-z0-9]{15,}',
 r'(?i)(?:token|password|secret)\s*[:=]\s*[\'"][A-Za-z0-9_-]{16,}',
 r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}']
for p in root.rglob('*'):
 if not p.is_file() or '.git' in p.parts:continue
 rel=p.relative_to(root)
 if not (str(rel) in top or (rel.parts[0] in {'R','docs','tests'} and p.suffix in exts)):bad.append(str(rel)+': not allowlisted')
 t=p.read_text(encoding='utf-8')
 if p.name=='check_distribution.py':continue # scanner patterns are not secrets
 if any(re.search(q,t) for q in patterns):bad.append(str(rel)+': sensitive-pattern match')
 if any(x in t.lower() for x in ['fjdctor','fjdoc','huadong','xiaoxin','jiuyi']):bad.append(str(rel)+': personal metadata')
assert not bad,'\n'.join(bad)
z=json.loads((root/'.zenodo.json').read_text())
assert z['upload_type']=='software' and z['license']=='mit'
print('PASS: code-only allowlist and sensitive-pattern checks')
