#!/usr/bin/env python3
"""Replay focus and 1000 q rows against separately enumerated frozen input."""
import collections,json,multiprocessing as mp,pathlib,sys,time
from driver import enumerate_stage,identity,limits
from replay_all_impl import verify_one
P=pathlib.Path
if __name__=='__main__':
 data=P(sys.argv[1]);b=P(sys.argv[2]);rs={};groups={}
 for f in sorted(b.glob('*.jsonl')):
  if not(f.stem.startswith('focus_') or f.stem=='q_fibre1000'):continue
  for l in f.open():
   r=json.loads(l);t=r.pop('task');r.pop('old_status',None);assert identity(t)==r['task_sha256'];assert r['id'] not in rs;rs[r['id']]=r;groups[r['id']]=f.stem
 found={};start=time.monotonic()
 def tasks():
  for s in ('B_std','B_E2','B_E3','q_other'):
   for t in enumerate_stage(data,s):
    if t['id'] in rs:found[t['id']]=1;yield t,rs[t['id']]
 ct=collections.defaultdict(collections.Counter)
 with mp.get_context('spawn').Pool(2,initializer=limits,initargs=(1024*1024**2,)) as pool:
  for i,status,why in pool.imap_unordered(verify_one,tasks(),chunksize=4):
   ct[groups[i]][status]+=1
   if status=='REJECT':print('REJECT',i,why,flush=True)
 assert set(found)==set(rs);print('FOCUS REPLAY',json.dumps(ct,sort_keys=True),'wall',time.monotonic()-start,'rows',len(rs),flush=True)
 if any(c.get('REJECT') for c in ct.values()):sys.exit(2)
