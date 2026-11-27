#!/usr/bin/env python3
"""Full qtable geometry and easy analytic bounds; does not trust stored LG labels."""
import json,collections,itertools,argparse,pathlib,time
from fractions import Fraction as Q
from label_exact import PI2,lower_quad,quad,ear_ranges

def convex_possible(lo,hi):
 _,s,_,_,vs=quad('std',list(lo)+[0,0],list(hi)+[1,1]);return all(b>0 for a,b in ear_ranges(vs)[:4])
def main():
 p=argparse.ArgumentParser();p.add_argument('data',type=pathlib.Path);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();st=time.monotonic();ct=collections.Counter();bad=[];partsct=collections.defaultdict(collections.Counter);roots=set()
 with open(a.out,'w',buffering=1) as f:
  for i,line in enumerate(open(a.data/'qtable/cells.jsonl'),1):
   r=json.loads(line);k=tuple(r['key']);assert k not in roots and len(k)==5 and k[0]==0;roots.add(k);ct['rows']+=1
   lo=[Q(z,16) for z in k[1:]];hi=[x+Q(1,16) for x in lo];L=Q(r['LB']);parts=r.get('lg',{}).get('parts',[])
   if parts:
    ct['lg_rows']+=1;sub=set();whole=False;loLs=[]
    for j,(l,h,Lp,met) in enumerate(parts):
     l=list(map(Q,l));h=list(map(Q,h));Lp=Q(Lp);loLs.append(Lp);partsct[met]['total']+=1
     if l==lo and h==hi:whole=True
     else:
      v=tuple((x-y)*32 for x,y in zip(l,lo))
      if any(x not in (0,1) for x in v) or any(b-a!=Q(1,32) for a,b in zip(l,h)) or v in sub:bad.append((i,j,'part geometry'))
      sub.add(v)
     # Independently prove a subset of all stored parts using analytic bounds.
     _,s,_,_,qs=quad('std',l+[0,0],h+[1,1]);rr=ear_ranges(qs);lb=lower_quad('std',s,qs,rr)
     if lb>=Lp or any(b<=0 for a,b in rr[:4]):partsct[met]['analytic_pass']+=1;f.write(json.dumps({'id':f'qstd/{i}/{j}','method':met,'status':'ANALYTIC_PROVED','L':str(Lp)})+'\n')
     else:partsct[met]['requires_more']+=1
    if L>min(loLs):bad.append((i,'root LB exceeds part'))
    if whole:
     if len(parts)!=1:bad.append((i,'whole and subparts'))
    else:
     for e in itertools.product((0,1),repeat=4):
      l=[x+Q(b,32) for x,b in zip(lo,e)];h=[x+Q(1,32) for x in l]
      if e not in sub and convex_possible(l,h):bad.append((i,'missing convex part',e))
   else:
    ct['closed_rows']+=1;_,s,_,_,qs=quad('std',lo+[0,0],hi+[1,1]);rr=ear_ranges(qs);lb=lower_quad('std',s,qs,rr)
    met=r['how'];partsct['root_'+met]['total']+=1
    if lb>=L or any(b<=0 for a,b in rr[:4]):partsct['root_'+met]['analytic_pass']+=1;f.write(json.dumps({'id':f'qstd/{i}','method':met,'status':'ANALYTIC_PROVED','L':str(L)})+'\n')
    else:partsct['root_'+met]['requires_more']+=1
   if i%10000==0:print('ROWS',i,'bad',len(bad),'wall',time.monotonic()-st,flush=True)
 rootdata=json.load(open(a.data/'std/roots.json'))['blocks'];expected={tuple(k):Q(L) for rows in rootdata.values() for k,*_,L in rows};print('ROOTS MATCH',roots==set(expected),len(roots),'BAD',len(bad),flush=True)
 print('COUNTS',json.dumps(ct),'PARTS',json.dumps(partsct),'wall',time.monotonic()-st,flush=True)
 print('BAD_LIST',json.dumps(bad),flush=True)
if __name__=='__main__':main()
