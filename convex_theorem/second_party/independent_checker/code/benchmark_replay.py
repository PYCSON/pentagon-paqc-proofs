#!/usr/bin/env python3
"""Convert actual benchmark rows to the normal replay format (sample only, missing explicit)."""
import argparse,collections,json,pathlib
from driver import STAGES,VERSION,ZIP_SHA,CHUNK,code_id,enumerate_stage,task_iter,chunks,identity,sha,canonical,atomic
P=pathlib.Path
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('data',type=P);p.add_argument('out',type=P);p.add_argument('bench',nargs='+',type=P);p.add_argument('--sample',type=int,default=200);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 full={s:sum(1 for _ in enumerate_stage(a.data,s)) for s in STAGES};totals={s:min(n,a.sample) for s,n in full.items()};summary={}
 atomic(a.out/'plan.json',dict(config=dict(version=VERSION,code=code_id(),zip_sha256=ZIP_SHA,sample=a.sample,engine='scipy',budget=64,chunk=CHUNK),totals=totals,full_totals=full))
 for s in STAGES:
  rr={}
  for d in a.bench:
   f=d/('bench_'+s+'.jsonl')
   if not f.exists():continue
   for line in f.open():
    r=json.loads(line);t=r.pop('task');r.pop('old_status',None);assert identity(t)==r['task_sha256'];rr[r['id']]=r
  ct=collections.Counter();seen=set()
  for i,b in enumerate(chunks(task_iter(a.data,s,a.sample))):
   path=a.out/s/f'{i:06d}.jsonl';rs=[]
   for t in b:
    if t['id'] in rr:
     r=rr[t['id']];assert r['task_sha256']==identity(t);rs.append(r);ct[r['status']]+=1;seen.add(r['id'])
   if rs:path.parent.mkdir(parents=True,exist_ok=True);path.write_text(''.join(canonical(r)+'\n' for r in rs))
  assert seen==set(rr);summary[s]=dict(total=totals[s],counts=dict(ct),missing=totals[s]-sum(ct.values()),note='BENCHMARK SAMPLE ONLY')
 atomic(a.out/'summary.json',summary)
 paths=sorted(p for p in a.out.rglob('*') if p.is_file());(a.out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(a.out)}\n' for p in paths));print(json.dumps(summary,sort_keys=True))
