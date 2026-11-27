#!/usr/bin/env python3
import argparse,collections,json,time
from fractions import Fraction as Q
from lg_exact import replay
from multiprocessing import Pool

def work(row):
 w=row.get('witness'); rid=row.get('id','unknown')
 if w is None:return rid,'NO_WITNESS'
 try:
  if w['chart']!=row['chart'] or w['j']!=row['j'] or list(map(Q,w['lo']))!=list(map(Q,row['lo'])) or list(map(Q,w['hi']))!=list(map(Q,row['hi'])) or Q(w['L'])<Q(row['L']):return rid,'BINDING_REJECT'
  return rid,'PASS' if replay(w) else 'REJECT'
 except Exception as e:return rid,'ERROR '+repr(e)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('files',nargs='+');p.add_argument('--jobs',type=int,default=1);p.add_argument('--sample',type=int,default=0);a=p.parse_args();st=time.monotonic();r=[json.loads(s) for f in a.files for s in open(f)];r=[x for x in r if x.get('witness')];r=r[:a.sample] if a.sample else r;c=collections.Counter()
 with Pool(a.jobs) as pool:
  for i,s in pool.imap_unordered(work,r,chunksize=1):c[s]+=1;print(i,s,flush=True)
 print('SUMMARY',dict(c),'wall',time.monotonic()-st,flush=True)
