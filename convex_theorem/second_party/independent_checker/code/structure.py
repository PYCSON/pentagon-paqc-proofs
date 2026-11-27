#!/usr/bin/env python3
"""Exact independent dyadic tiling + root-domain coverage. No producer imports."""
import argparse,collections,hashlib,itertools,json,pathlib,time
from fractions import Fraction as Q
from multiprocessing import Pool
from math import prod
S=1<<40
def toint(v):
 a,b=Q(v).as_integer_ratio();assert S%b==0;return a*(S//b)
def rect(lo,hi):return tuple(map(toint,lo)),tuple(map(toint,hi))
def tile(root,bs):
 lo,hi=root
 if not bs:return False
 for a,b in bs:
  if not all(l<=x<y<=h for l,x,y,h in zip(lo,a,b,hi)):return False
 if sum(prod(y-x for x,y in zip(a,b)) for a,b in bs)!=prod(y-x for x,y in zip(lo,hi)):return False
 def rec(lo,hi,bs):
  if len(bs)==1:return bs[0]==(lo,hi)
  for i in range(len(lo)):
   if (lo[i]+hi[i])%2:continue
   m=(lo[i]+hi[i])//2;left=[];right=[]
   for a,b in bs:
    if b[i]<=m:left.append((a,b))
    elif a[i]>=m:right.append((a,b))
    else:break
   else:
    if left and right:return rec(lo,hi[:i]+(m,)+hi[i+1:],left) and rec(lo[:i]+(m,)+lo[i+1:],hi,right)
  return False
 return rec(lo,hi,bs)
def check_block(task):
 chart,p,roots=task;groups=collections.defaultdict(list);counts=collections.Counter();errors=[];d=json.load(open(p/'done.json'))
 for name,sha in d['sha256'].items():
  b=(p/name).read_bytes() if (p/name).exists() else b'' if name=='corners.jsonl' else None
  if b is None or hashlib.sha256(b).hexdigest()!=sha:errors.append('hash '+name)
 for name in ('leaves.jsonl','cert.jsonl'):
  for line in open(p/name):
   r=json.loads(line)
   if name=='leaves.jsonl':k,j,lo,hi,label=r
   else:k,j,lo,hi,label=r['key'],r['j'],r['lo'],r['hi'],'C';assert r['how']=='LG' and r['status']=='certified'
   groups[tuple(k),j].append(rect(lo,hi));counts[label]+=1
 expected={(tuple(k),j):rect(lo+[0,0],hi+[1,1]) for k,lo,hi,lb in roots for j in [0,3]}
 if set(groups)!=set(expected):errors.append('root keys')
 for key,root in expected.items():
  if not tile(root,groups.get(key,[])):errors.append('tile '+str(key))
 return chart,p.name,len(expected),dict(counts),errors
C=list(map(Q,[1.3090169943749475,0.9510565162951536,.5,1.5388417685876268,1-1.3090169943749475,0.9510565162951536]))
def rad(lo,hi):
 a=[Q(x)-c for x,c in zip(lo,C)];b=[Q(x)-c for x,c in zip(hi,C)];return sum(min(x*x,y*y) if x*y>0 else 0 for x,y in zip(a,b)),sum(max(x*x,y*y) for x,y in zip(a,b))
def skip(lo,hi):
 x,y=rad(lo,hi)
 if Q(hi[2])<=Q(1,2):return 'mirror'
 if x>=Q(9,25):return 'far'
 if y<Q(1,400):return 'inner'
 return 'LG'
def near(path,out):
 root=[Q(round(c*1024),1024)-Q(5,8) for c in C];stored=set();ct=collections.Counter();errors=[];inner=0
 for p in sorted((path/'blocks').glob('*.jsonl')):
  f=open(p);meta=json.loads(next(f));lo,hi=meta['lo'],meta['hi'];k=tuple(map(int,p.stem));a=[r+i*Q(5,16) for r,i in zip(root,k)];b=[x+Q(5,16) for x in a]
  if list(map(Q,lo))!=a or list(map(Q,hi))!=b:errors.append(p.stem+' root mismatch')
  stored.add(k);bs=[]
  for line in f:
   l,h,lab,N=json.loads(line);ct[lab]+=1;bs.append(rect(l,h));s=skip(l,h)
   if lab!=s:errors.append(p.stem+' label '+lab+' != '+s)
   if lab=='LG':
    x,y=rad(l,h);inner+=x<=Q(1,400)<=y
  if not tile(rect(lo,hi),bs):errors.append(p.stem+' tiling')
 wanted=set()
 for k in itertools.product(range(4),repeat=6):
  a=[r+i*Q(5,16) for r,i in zip(root,k)];b=[x+Q(5,16) for x in a]
  if skip(a,b)=='LG':wanted.add(k)
 print('NEAR',len(stored),len(wanted),'missing',len(wanted-stored),'extra',len(stored-wanted),'labels',ct,'touch',inner,'errors',len(errors),flush=True);out.write(json.dumps({'near':dict(ct),'blocks':len(stored),'wanted':len(wanted),'missing':[list(k) for k in sorted(wanted-stored)],'extra':[list(k) for k in sorted(stored-wanted)],'touching':inner,'errors':errors})+'\n')
def qverts(chart,w):
 a,b,c,d=w
 if chart=='std':return (16*a,16*b),(16*c,16*d)
 if chart=='E2':return (256+c*(a-16),c*b+16*d),(16*a,16*b)
 return (16*a,16*b),(c*a,c*b+16*d)
def root_check(task):
 chart,stored=task;rg={'std':(range(32),range(16),range(-16,16),range(16)),'E2':(range(-16,16),range(4,16),range(-72,72),range(2)),'E3':(range(4,32),range(16),range(-64,64),range(2))}[chart];missing=[];extra=[];total=0;need=0
 for k in itertools.product(*rg):
  total+=1
  vs=[qverts(chart,tuple(x+e for x,e in zip(k,es))) for es in itertools.product((0,1),repeat=4)]
  vals=[(b,d,(a-256)*(d-b)-b*(c-a),a*d-b*c) for (a,b),(c,d) in vs]
  required=all(max(col)>0 for col in zip(*vals))
  if required and chart=='std':
   for edge in ([(a-256,b) for (a,b),(c,d) in vs],[(c-a,d-b) for (a,b),(c,d) in vs],[(c,d) for (a,b),(c,d) in vs]):
    lb=sum(0 if min(col)<=0<=max(col) else min(abs(min(col)),abs(max(col)))**2 for col in zip(*edge))
    if lb>256**2:required=False;break
  need+=required
  if required and k not in stored:missing.append(k)
  if not required and k in stored:extra.append(k)
 return dict(chart=chart,total=total,needed=need,stored=len(stored),missing=missing,extra=extra)
def main():
 p=argparse.ArgumentParser();p.add_argument('far',type=pathlib.Path);p.add_argument('near',type=pathlib.Path);p.add_argument('--jobs',type=int,default=4);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();t=time.monotonic();tasks=[];roots=[]
 for chart in ('std','E2','E3'):
  data=json.load(open(a.far/chart/'roots.json'))['blocks'];stored=set()
  assert set(data)=={x.name for x in (a.far/chart/'blocks').iterdir() if x.is_dir()}
  for name,rows in data.items():
   for key,lo,hi,L in rows:
    k=tuple(key[1:] if chart=='std' else key);assert k not in stored;stored.add(k);assert list(map(Q,lo))==[Q(x,16) for x in k] and list(map(Q,hi))==[Q(x+1,16) for x in k]
   tasks.append((chart,a.far/chart/'blocks'/name,rows))
  roots.append((chart,stored))
 counts=collections.Counter();errs=0;columns=0
 with open(a.out,'w',buffering=1) as out,Pool(a.jobs) as pool:
  for i,(chart,name,n,ct,e) in enumerate(pool.imap_unordered(check_block,tasks,chunksize=2),1):
   columns+=n;counts.update(ct);errs+=len(e);out.write(json.dumps(dict(chart=chart,block=name,roots=n,counts=ct,errors=e))+'\n')
   if i%100==0:print('blocks',i,'columns',columns,'boxes',sum(counts.values()),'errors',errs,flush=True)
  print('FAR',len(tasks),columns,sum(counts.values()),dict(counts),'errors',errs,'seconds',time.monotonic()-t,flush=True)
  for r in pool.imap_unordered(root_check,roots):
   print('ROOT',r['chart'],r['total'],r['needed'],r['stored'],'missing',len(r['missing']),'extra',len(r['extra']),flush=True);out.write(json.dumps(r)+'\n')
  near(a.near,out)
 print('WALL',time.monotonic()-t,flush=True)
if __name__=='__main__':main()
