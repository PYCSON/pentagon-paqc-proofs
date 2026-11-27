#!/usr/bin/env python3
"""C not explicitly supplied as a binary64 list: verify robust to its stated chart uncertainty."""
from fractions import Fraction as Q
from math import isqrt
from itertools import product
from pathlib import Path
import json,time,collections
from structure import C
B=1<<100
def sqiv(x):
 k=isqrt((x*B*B).__floor__());return Q(k,B),Q(k+1,B)
s5=sqiv(Q(5));x0=((3+s5[0])/4,(3+s5[1])/4);y1=(sqiv((5+s5[0])/8)[0],sqiv((5+s5[1])/8)[1]);y2=(sqiv((5-s5[1])/8)[0],sqiv((5-s5[0])/8)[1]);top=(y1[0]+y2[0],y1[1]+y2[1]);vals=[x0,y1,(Q(1,2),Q(1,2)),top,(1-x0[1],1-x0[0]),y1]
assert all(max(abs(c-l),abs(c-h))<=Q(1,1<<50) for c,(l,h) in zip(C,vals));eps=Q(1,1<<49)
# |C_archived - C_here| <= eps, exact regular-chart radicals checked above.
def radii(lo,hi):
 a=[Q(x)-c for x,c in zip(lo,C)];b=[Q(x)-c for x,c in zip(hi,C)];mn=[min(abs(x),abs(y)) if x*y>0 else Q(0) for x,y in zip(a,b)];mx=[max(abs(x),abs(y)) for x,y in zip(a,b)];return sum(max(0,x-eps)**2 for x in mn),sum((x+eps)**2 for x in mn),sum(max(0,x-eps)**2 for x in mx),sum((x+eps)**2 for x in mx)
def label(lo,hi):
 if Q(hi[2])<=Q(1,2):return 'mirror'
 l,L,h,H=radii(lo,hi)
 if l>=Q(9,25):return 'far'
 if H<Q(1,400):return 'inner'
 if L<Q(9,25) and h>=Q(1,400):return 'LG'
 return 'UNCERTAIN'
st=time.monotonic();p=Path(__import__('sys').argv[1] if len(__import__('sys').argv)>1 else '/mnt/data/p24audit/data/near')/'blocks';ct=collections.Counter();bad=[];inner=0;uncertain=0;stored=set()
for f in sorted(p.glob('*.jsonl')):
 lines=open(f);meta=json.loads(next(lines));stored.add(tuple(map(int,f.stem)))
 for i,line in enumerate(lines,1):
  lo,hi,lab,N=json.loads(line);r=label(lo,hi);ct[lab]+=1
  if r!=lab:bad.append((f.stem,i,lab,r))
  if lab=='LG':
   l,L,h,H=radii(lo,hi);inner+=L<=Q(1,400)<=h;uncertain+=(l<=Q(1,400)<=H) and not (L<=Q(1,400)<=h)
root=[Q(round(c*1024),1024)-Q(5,8) for c in C];assert all(abs(c-(a+Q(5,8)))+eps<Q(1,2048) for a,c in zip(root,C));wanted=set()
for k in product(range(4),repeat=6):
 lo=[a+i*Q(5,16) for a,i in zip(root,k)];hi=[x+Q(5,16) for x in lo];r=label(lo,hi)
 if r=='LG':wanted.add(k)
 elif r=='UNCERTAIN':bad.append(('root',k,r))
print('ROBUST C',list(map(str,C)),'eps',str(eps),'leaves',dict(ct),'roots',len(wanted),'missing',len(wanted-stored),'extra',len(stored-wanted),'touch',inner,'touch_uncertain',uncertain,'errors',len(bad),'wall',time.monotonic()-st);print('bad',bad[:30]);assert not bad and wanted==stored
