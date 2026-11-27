"""Independent exact predicate and bounded near-only supplementary refit (r3 untouched)."""
from fractions import Fraction as Q
import time
from proofs24b import mid_split,refit_lg
import driver as D
HEX=('0x1.4f1bbcdcbfa53p+0','0x1.e6f0e134454ffp-1','0x1.0p-1','0x1.89f188bdcd7aep+0','-0x1.3c6ef372fe94cp-2','0x1.e6f0e134454ffp-1')
C=tuple(Q(float.fromhex(x)) for x in HEX)
OUT2=(Q(3,5)+Q(1,2**40))**2
IDS=frozenset('near/'+x for x in '023112/23 023121/14 132012/151 132012/26 132012/47 132012/7 132022/157 132022/72 132102/10 133022/30 133022/40 232012/19 233021/2'.split())
def inf_d2(lo,hi):
 if len(lo)!=6 or len(hi)!=6:raise ValueError('near dimension')
 s=Q(0)
 for l,h,c in zip(map(Q,lo),map(Q,hi),C):
  if l>h:raise ValueError('reversed box')
  s+=max(l-c,c-h,Q(0))**2
 return s
def outside(lo,hi):return inf_d2(lo,hi)>=OUT2

def search(t,engine='scipy',budget=256):
 if t['id'] not in IDS or t['stage']!='near' or t['chart']!='near':raise ValueError('near13 only')
 start=time.monotonic();n=min(12,t['N']+2);r={'id':t['id'],'task_sha256':D.identity(t),'status':'UNPROVED'};lg=outs=attempts=0;reason='depth/budget exhausted';
 def rec(lo,hi,dep):
  nonlocal lg,outs,attempts,reason
  if outside(lo,hi):outs+=1;return {'outside':True}
  if lg>=budget or dep>14:return None
  # Refit on the actual retained node, using one exact LG leaf only. Do not spend an LG
  # trial at a huge node before six axes have been cut once. All children first test outside.
  if dep>=6:
   attempts+=1;u=dict(t,lo=lo,hi=hi,N=n)
   w,info=refit_lg(u,engine,1)
   if w is not None:lg+=1;return w
   reason=info.get('reason','LG positivity not established')
  if dep==14:return None
  # longest relative side: first six cuts divide every original axis once
  axis=max(range(6),key=lambda i:(Q(hi[i])-Q(lo[i]))/(Q(t['hi'][i])-Q(t['lo'][i])))
  a,b=mid_split(lo,hi,axis);left=rec(*a,dep+1)
  if left is None:return None
  right=rec(*b,dep+1)
  if right is None:return None
  return {'axis':axis,'left':left,'right':right}
 try:
  w=rec(list(map(Q,t['lo'])),list(map(Q,t['hi'])),0)
  if w is not None and lg>0:r.update(status='PROVED',witness=w,reason='near-only coverage')
  else:r['reason']=reason
 except Exception as ex:r.update(status='ERROR',reason=repr(ex)[:400])
 r.update(parts=lg,outside_leaves=outs,attempts=attempts,seconds=time.monotonic()-start)
 if len(D.canonical(r).encode())+1>D.MAX_ROW:r.pop('witness',None);r.update(status='UNPROVED',reason='row cap exceeded')
 return r
