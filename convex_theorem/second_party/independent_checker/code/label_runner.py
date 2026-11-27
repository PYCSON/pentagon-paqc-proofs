#!/usr/bin/env python3
import argparse,collections,json,pathlib,time,hashlib,random,itertools
from multiprocessing import Pool,TimeoutError
from fractions import Fraction as Q
from label_exact import check

def split_all(lo,hi):
 m=[(Q(a)+Q(b))/2 for a,b in zip(lo,hi)]
 for es in itertools.product((0,1),repeat=6):yield [a if not e else x for a,x,e in zip(lo,m,es)],[x if not e else b for x,b,e in zip(m,hi,es)]
def recheck(chart,j,lo,hi,lab,table,depth):
 r=check(chart,j,lo,hi,lab,table)
 if r:return r,False
 if not depth:return None,False
 conditional=False
 for a,b in split_all(lo,hi):
  r,_=recheck(chart,j,a,b,lab,table,depth-1)
  if not r:return None,False
  conditional|=r=='CONDITIONAL_TABLE'
 return ('CONDITIONAL_TABLE' if conditional else 'PROVED'),True

def work(t):
 chart,block,path,roots,labels,sample,depth=t;ct=collections.defaultdict(collections.Counter);bad=[];start=time.monotonic();rows=[json.loads(l) for l in open(path)]
 if sample:
  by=collections.defaultdict(list)
  for i,r in enumerate(rows,1):
   if r[4] in labels:by[r[4]].append((i,r))
  seq=[]
  for lab,v in by.items():random.Random(20261122).shuffle(v);seq+=v[:sample]
 else:seq=[(i,r) for i,r in enumerate(rows,1) if r[4] in labels]
 for i,(k,j,lo,hi,lab) in seq:
  ct[lab]['checked']+=1
  try:r,split=recheck(chart,j,lo,hi,lab,roots[tuple(k)],depth);reason='criterion did not close'
  except Exception as e:r=None;split=False;reason=repr(e)
  if r:
   ct[lab]['split' if split else 'direct']+=1
   if r=='CONDITIONAL_TABLE':ct[lab]['conditional']+=1
  else:ct[lab]['unproved']+=1;bad.append(dict(chart=chart,block=block,line=i,key=k,j=j,lo=lo,hi=hi,label=lab,reason=reason))
 return dict(chart=chart,block=block,counts=ct,unproved=bad,seconds=time.monotonic()-start)

def main():
 p=argparse.ArgumentParser();p.add_argument('data',type=pathlib.Path);p.add_argument('--jobs',type=int,default=1);p.add_argument('--blocks',nargs='*');p.add_argument('--sample',type=int,default=0,help='per label per block; 0 means all');p.add_argument('--labels',nargs='+',default=['V','E2','E3','X','N','S','R','T','TX','O']);p.add_argument('--depth',type=int,default=2,choices=[0,1,2]);p.add_argument('--seconds',type=int,default=0);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();assert 1<=a.jobs<=4
 jobs=[]
 for chart in ('std','E2','E3'):
  roots=json.load(open(a.data/chart/'roots.json'))['blocks']
  for block,rs in sorted(roots.items()):
   if a.blocks and block not in a.blocks and chart+'/'+block not in a.blocks:continue
   p=a.data/chart/'blocks'/block/'leaves.jsonl';r={tuple(k):L for k,lo,hi,L in rs};jobs.append((chart,block,p,r,set(a.labels),a.sample,a.depth))
 config=hashlib.sha256(repr((sorted(a.labels),a.sample,a.depth)).encode()).hexdigest();done={}
 if a.out.exists():
  for l in open(a.out):
   r=json.loads(l);assert r['config']==config;done[r['chart'],r['block']]=r
 todo=[t for t in jobs if (t[0],t[1]) not in done];print('BLOCKS',len(jobs),'PENDING',len(todo),'LABELS',a.labels,flush=True)
 st=time.monotonic()
 with open(a.out,'a',buffering=1) as f,Pool(a.jobs) as pool:
  it=pool.imap_unordered(work,todo,chunksize=1)
  for i in range(len(todo)):
   try:r=it.next(timeout=max(.01,a.seconds-(time.monotonic()-st)) if a.seconds else None)
   except TimeoutError:print('TIME_LIMIT',flush=True);break
   r['config']=config;done[r['chart'],r['block']]=r;f.write(json.dumps(r,separators=(',',':'))+'\n')
   if i%25==0:print('BLOCK',i+1,r['chart'],r['block'],r['counts'],flush=True)
   if a.seconds and time.monotonic()-st>=a.seconds:break
 acc=collections.defaultdict(collections.Counter)
 for (chart,block),r in done.items():
  for lab,c in r['counts'].items():acc[chart+'/'+lab].update(c)
 print('SUMMARY',json.dumps(acc,sort_keys=True),'completed_blocks',len(done),'wall',time.monotonic()-st,flush=True)
if __name__=='__main__':main()
