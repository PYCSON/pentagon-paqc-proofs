#!/usr/bin/env python3
"""Independent candidate generator + strict verifier. Resume by input identity and config.
Search uses NumPy/SciPy; proof replay lg_exact.py uses stdlib only."""
import argparse, pathlib, json, time, random, hashlib, os, collections
from multiprocessing import Pool
from fractions import Fraction as Q
from lg_exact import TARGET
C=list(map(Q,[1.3090169943749475,0.9510565162951536,.5,1.5388417685876268,1-1.3090169943749475,0.9510565162951536]))
def bounds(lo,hi):
 a=[Q(x)-c for x,c in zip(lo,C)]; b=[Q(x)-c for x,c in zip(hi,C)];return sum(min(x*x,y*y) if x*y>0 else 0 for x,y in zip(a,b)),sum(max(x*x,y*y) for x,y in zip(a,b))
def near_tasks(path):
 for p in sorted((path/'blocks').glob('*.jsonl')):
  with open(p) as f:
   next(f)
   for i,line in enumerate(f,1):
    lo,hi,lab,N=json.loads(line)
    if lab=='LG':yield dict(id=f'near/{p.stem}/{i}',chart='near',lo=lo,hi=hi,j=3,N=N,L=str(TARGET),block=p.stem)
def far_tasks(path):
 for chart in ('std','E2','E3'):
  for p in sorted((path/chart/'blocks').glob('*/cert.jsonl')):
   for i,line in enumerate(open(p),1):
    q=json.loads(line);yield dict(id=f'{chart}/{p.parent.name}/{i}',chart=chart,lo=q['lo'],hi=q['hi'],j=q['j'],N=q['N'],L=str(TARGET),block=p.parent.name,key=q['key'])
def q_tasks(path):
 for i,line in enumerate(open(path/'qtable/cells.jsonl'),1):
  r=json.loads(line)
  for j,(lo,hi,L,met) in enumerate(r.get('lg',{}).get('parts',[])):
   yield dict(id=f'qstd/{i}/{j}',chart='qstd',lo=lo,hi=hi,j=3,N=r['lg']['N'],L=str(Q(L)),block=str(i),method=met,key=r['key'])
def work(t):
 st=time.monotonic(); row=dict(t)
 try:
  if t.get('engine')=='stdlib':
   from lg_stdlib import certify
  else:
   from lg_search import certify
  L=Q(t['L']);w,info=certify(t['chart'],t['lo'],t['hi'],t['j'],t['N'],L,budget=t['budget']);row.update(info);row['status']='PROVED' if w is not None else 'UNPROVED';row['witness']=w
 except Exception as e:row.update(status='ERROR',reason=repr(e))
 row['seconds']=time.monotonic()-st;return row

def main():
 p=argparse.ArgumentParser();p.add_argument('data',type=pathlib.Path);p.add_argument('--kind',choices=['near','far','q'],required=True);p.add_argument('--jobs',type=int,default=1);p.add_argument('--blocks',nargs='*');p.add_argument('--charts',nargs='*');p.add_argument('--sample',type=int,default=0);p.add_argument('--priority',choices=['inner','largest','nonconvex','perblock','random','first','all'],default='all');p.add_argument('--engine',choices=['scipy','stdlib'],default='scipy');p.add_argument('--per-block',type=int,default=3);p.add_argument('--retry-unproved',action='store_true');p.add_argument('--budget',type=int,default=64);p.add_argument('--out',type=pathlib.Path,required=True);p.add_argument('--seconds',type=int,default=0);a=p.parse_args();assert 1<=a.jobs<=4 and 1<=a.budget<=64
 tasks=list({'near':near_tasks,'far':far_tasks,'q':q_tasks}[a.kind](a.data)); print('ENUMERATED',len(tasks),flush=True)
 if a.blocks:tasks=[x for x in tasks if x['block'] in a.blocks or x['chart']+'/'+x['block'] in a.blocks]
 if a.charts:tasks=[x for x in tasks if x['chart'] in a.charts]
 if a.kind=='q':tasks=[x for x in tasks if x['method']=='lg']
 if a.priority=='inner':tasks=[x for x in tasks if (lambda b:b[0]<=Q(1,400)<=b[1])(bounds(x['lo'],x['hi']))];print('INNER',len(tasks),flush=True)
 if a.priority=='largest':
  by={}
  from math import prod
  for t in tasks:
   k=(t['chart'],t['block']);v=prod(Q(y)-Q(x) for x,y in zip(t['lo'],t['hi']))
   if k not in by or v>by[k][0]:by[k]=(v,t)
  tasks=[t for v,t in by.values()]
 if a.priority=='perblock':
  by=collections.defaultdict(list)
  for t in tasks:by[t['chart'],t['block']].append(t)
  tasks=[]
  for v in by.values():random.Random(20261122).shuffle(v);tasks+=v[:a.per_block]
 if a.priority=='nonconvex':
  if a.kind!='near':raise ValueError('nonconvex selection is near only')
  from itertools import product
  from label_exact import cross,vecsub
  random.Random(20261122).shuffle(tasks);selected=[]
  for t in tasks:
   lo,hi=t['lo'],t['hi'];
   for z in product(*zip(map(Q,lo),map(Q,hi))):
    V=[(0,0),(1,0),*[(z[i],z[i+1]) for i in (0,2,4)]]
    if any(cross(vecsub(V[i],V[i-1]),vecsub(V[(i+1)%5],V[i]))<0 for i in range(5)):
     selected.append(t);break
   if a.sample and len(selected)>=a.sample:break
  tasks=selected
 if a.priority=='random':random.Random(20261122).shuffle(tasks)
 if a.sample:tasks=tasks[:a.sample]
 for t in tasks:
  t['budget']=a.budget
  if a.engine!='scipy':t['engine']=a.engine
  t['input_id']=hashlib.sha256(json.dumps(t,sort_keys=True).encode()).hexdigest()
 a.out.parent.mkdir(parents=True,exist_ok=True);done={}
 if a.out.exists():
  for line in open(a.out):
   try:r=json.loads(line);done[r['input_id']]=r
   except Exception:raise RuntimeError('truncated result line; repair before resume')
 pending=[t for t in tasks if t['input_id'] not in done or (a.retry_unproved and done[t['input_id']]['status']!='PROVED')];print('SELECTED',len(tasks),'PENDING',len(pending),flush=True); st=time.monotonic()
 with open(a.out,'a',buffering=1) as f,Pool(a.jobs) as pool:
  it=pool.imap_unordered(work,pending,chunksize=1)
  for _ in range(len(pending)):
   try:r=it.next(timeout=max(.01,a.seconds-(time.monotonic()-st)) if a.seconds else None)
   except TimeoutError:break
   except __import__('multiprocessing').TimeoutError:print('TIME_LIMIT',flush=True);break
   f.write(json.dumps(r,separators=(',',':'))+'\n');done[r['input_id']]=r; print(r['id'],r['status'],r.get('parts'),r.get('lower_center'),round(r['seconds'],4),r.get('reason'),flush=True)
   if a.seconds and time.monotonic()-st>=a.seconds:break
 print('SUMMARY',collections.Counter(done[t['input_id']]['status'] if t['input_id'] in done else 'NOT_RUN' for t in tasks),'wall',time.monotonic()-st,flush=True)
if __name__=='__main__':main()
