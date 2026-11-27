#!/usr/bin/env python3
"""Independent full qtable non-LG proof runner. Uses analytic bounds or strict PS fibre test.
Does not count archived `how`/`closed` labels as proof. Exact rational checks only."""
import argparse,json,collections,pathlib,time,random,hashlib
from fractions import Fraction as Q
from multiprocessing import Pool,TimeoutError
from label_exact import quad,ear_ranges,lower_quad,fibre

def tasks(P):
 for i,line in enumerate(open(P/'qtable/cells.jsonl'),1):
  r=json.loads(line);k=r['key'];L=Q(r['LB']);parts=r.get('lg',{}).get('parts',[])
  if parts:
   for j,(a,b,l,m) in enumerate(parts):
    if m!='lg':yield dict(id=f'qstd/{i}/{j}',lo=a,hi=b,L=str(Q(l)),method=m,block=str(i),key=k)
  else:
   a=[Q(x,16) for x in k[1:]];b=[x+Q(1,16) for x in a];yield dict(id=f'qstd/{i}',lo=list(map(str,a)),hi=list(map(str,b)),L=str(L),method='root_'+r['how'],block=str(i),key=k)
def work(r):
 st=time.monotonic();x=dict(r)
 try:
  _,s,_,_,qs=quad('std',list(map(Q,r['lo']))+[0,0],list(map(Q,r['hi']))+[1,1]);rr=ear_ranges(qs);L=Q(r['L']);lb=lower_quad('std',s,qs,rr)
  if any(b<=0 for a,b in rr[:4]) or lb>=L:x['status']='ANALYTIC_PROVED'
  else:
   vs=[((0,0),(s,0),a,b) for a,b,_,_ in qs];areas=[q[3] for q in qs];x['status']='FIBRE_PROVED' if fibre(vs,areas,L) else 'UNPROVED'
 except Exception as e:x.update(status='ERROR',reason=repr(e))
 x['seconds']=time.monotonic()-st;return x
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('data',type=pathlib.Path);p.add_argument('--jobs',type=int,default=1);p.add_argument('--sample',type=int,default=0);p.add_argument('--blocks',nargs='*');p.add_argument('--seconds',type=int,default=0);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();assert 1<=a.jobs<=4
 ts=list(tasks(a.data));ts=[t for t in ts if not a.blocks or t['block'] in a.blocks];random.Random(20261122).shuffle(ts);ts=ts[:a.sample] if a.sample else ts;done={}
 if a.out.exists():
  for l in open(a.out):r=json.loads(l);done[r['input_id']]=r
 for t in ts:t['input_id']=hashlib.sha256(json.dumps(t,sort_keys=True).encode()).hexdigest()
 pending=[t for t in ts if t['input_id'] not in done];st=time.monotonic();print('SELECTED',len(ts),'PENDING',len(pending),flush=True)
 with open(a.out,'a',buffering=1) as f,Pool(a.jobs) as pool:
  it=pool.imap_unordered(work,pending,chunksize=1)
  for i in range(len(pending)):
   try:r=it.next(timeout=max(.01,a.seconds-(time.monotonic()-st)) if a.seconds else None)
   except TimeoutError:print('TIME_LIMIT',flush=True);break
   f.write(json.dumps(r)+'\n');done[r['input_id']]=r;print(r['id'],r['status'],round(r['seconds'],4),flush=True)
   if a.seconds and time.monotonic()-st>=a.seconds:break
 print('SUMMARY',collections.Counter(done[t['input_id']]['status'] if t['input_id'] in done else 'NOT_RUN' for t in ts),'wall',time.monotonic()-st,flush=True)
