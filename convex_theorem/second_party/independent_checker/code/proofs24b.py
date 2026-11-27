"""Independent PENT24b proof generation/replay. No producer code. Numerical LG proposals only.
All accepted inequalities are exact rational/integer calculations. Label witnesses record the proof route.
"""
from fractions import Fraction as Q
from itertools import product
from functools import lru_cache
import time
import label_exact as E
import lg_exact as G

def mid_split(lo,hi,i):
 lo=list(map(Q,lo));hi=list(map(Q,hi));m=(lo[i]+hi[i])/2;a=hi.copy();a[i]=m;b=lo.copy();b[i]=m
 return (lo,a),(b,hi)

def refit_lg(t,engine='scipy',budget=64):
 # Each child solves its own first-order family at its own centre and scale.
 if engine=='scipy':import lg_search as S
 else:import lg_stdlib as S
 chart,j,n,L=t['chart'],t['j'],t['N'],Q(t['L']); attempts=0; last={}; maxdepth=10
 def rec(lo,hi,b,depth):
  nonlocal attempts,last
  attempts+=1
  if attempts>2*budget-1:return None,0
  try:
   u,chi=S.proposal(chart,lo,hi,j,n,order=1)
   ds,pd,pc,lb,ub,A,B=G.assemble(chart,lo,hi,j,n,u,chi,L,with_base=True);ps=ds+[pd,pc]
   last={'lower_center':float(lb) if lb is not None else None,'ritz_center':float(ub) if ub is not None else None,'terms_PD':len(pd),'terms_PC':len(pc)}
   d=len(lo);target=L
   def cand(ps):
    ms=[];bad=[]
    for i,p in enumerate(ps):
     m,z=S.choose_mu(p,d);ms.append(m)
     if z is None:bad.append(i)
    return ms,bad
   ms,bad=cand(ps)
   if bad and chart.startswith('q') and all(i==len(ps)-1 for i in bad):
    # Nonmonotone sufficient test: search directly for L' >= stored target.
    for k in (1,2,4,8,16,32,63):
     L1=L*(1+Q(k,16384));pp=G.sub(G.scale(A,L1.denominator),G.scale(B,L1.numerator));m,z=S.choose_mu(pp,d)
     if z is not None:target=L1;ps=ds+[pd,pp];ms[-1]=m;bad=[];break
   if not bad:
    w={'chart':chart,'lo':list(map(str,map(Q,lo))),'hi':list(map(str,map(Q,hi))),'j':j,'N':n,'L':str(target),'u':u,'chi':chi,'digest':G.digest(ps),'tree':{'mu':ms}}
    if not G.replay(w):raise AssertionError('LG kernel self-replay rejected')
    return {'lg':w},1
   if b<2 or depth>=maxdepth:return None,0
   # Only heuristic: score based on failed polynomial variation. Acceptance never uses it.
   scores=[sum(sum(abs(v)/max(1,abs(ps[k].get(0,0))) for e,v in ps[k].items() if (e>>(5*i))&31) for k in bad) for i in range(d)]
   axis=max(range(d),key=lambda i:scores[i])
   # Do not retain the parent polynomial dictionaries while refitting children.
   del ps,ds,pd,pc,A,B,u,chi
   if "pp" in locals():del pp
  except (ValueError,ArithmeticError) as ex:
   last={'reason':type(ex).__name__+': '+str(ex)[:300]}
   if b<2 or depth>=maxdepth:return None,0
   axis=max(range(len(lo)),key=lambda i:Q(hi[i])-Q(lo[i]))
  (a,c),(d,e)=mid_split(lo,hi,axis)
  left,nl=rec(a,c,b-1,depth+1)
  if left is None:return None,0
  right,nr=rec(d,e,b-nl,depth+1)
  if right is None:return None,0
  return {'axis':axis,'left':left,'right':right},nl+nr
 w,nparts=rec(t['lo'],t['hi'],budget,0);last.update(parts=nparts,attempts=attempts)
 return w,last

def lg_bound(t,w):
 def walk(w,lo,hi,depth=0):
  if depth>12:return 10000
  if 'lg' in w:
   if set(w)!={'lg'}:return 10000
   p=w['lg']
   if 'vertices' in p or p['chart']!=t['chart'] or p['j']!=t['j'] or p['N']!=t['N'] or Q(p['L'])<Q(t['L']) or list(map(Q,p['lo']))!=list(map(Q,lo)) or list(map(Q,p['hi']))!=list(map(Q,hi)):return 10000
   return 1 if G.replay(p) else 10000
  if set(w)!={'axis','left','right'} or type(w['axis']) is not int or not 0<=w['axis']<len(lo):return 10000
  a,b=mid_split(lo,hi,w['axis']);return walk(w['left'],*a,depth+1)+walk(w['right'],*b,depth+1)
 return walk(w,t['lo'],t['hi'])<=64

def rect_strip(vs,areas,target):
 # Concavity in each coordinate: pi^2 n^2 area2 - 2 target (projection difference)^2.
 for n in E.dirs(vs):
  n2=E.norm2(n)
  if n2 and all(E.PI2*n2*a>=2*target*(max(pr)-min(pr))**2 for a,v in zip(areas,vs) for pr in [[E.dot(n,x) for x in v]]) and any(a>0 for a in areas):
   # For strict inequality at every member use strict at positive-area corners;
   # if zero-area corners also zero on both sides, use strict margin via boosted target instead.
   return True
 return False

# More accurate step-minorant PS solution, with EXACT hull-breakpoint segments.
I=1<<80

def enc(x):
 x=Q(x)*I;return x.__floor__(),x.__ceil__()
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
 l,h=E.rng(x*y for x in a for y in b);return l//I,-((-h)//I)
def div(a,k):return a[0]//k,-((-a[1])//k)
def trans(z):
 if abs(z)>4:return None
 zz=enc(z);one=(I,I);C=S=tc=ts=one
 for k in range(1,21):
  tc=div(mul(tc,zz),(2*k-1)*(2*k));ts=div(mul(ts,zz),(2*k)*(2*k+1));C=add(C,tc);S=add(S,ts)
 # Both tails < 2*4^21/42! < 2^-80; directed rounding above.
 return (C[0]-1,C[1]+1),(S[0]-1,S[1]+1)

def hull_setup(H,d):
 n=E.norm2(d)
 if not n:return None
 xx=[E.dot(p,d) for p in H];yy=[E.cross(d,p) for p in H];w=max(xx)-min(xx);h=max(yy)-min(yy)
 if not w or not h:return None
 vv=[(Q(x-min(xx),w),Q(y-min(yy),w)) for x,y in zip(xx,yy)]
 U=[];L=[]
 for a,b in zip(vv,vv[1:]+vv[:1]):
  if a[0]==b[0]:continue
  r=(b[1]-a[1])/(b[0]-a[0]);(U if b[0]<a[0] else L).append((r,a[1]-r*a[0]))
 if not U or not L:return None
 return n,w,h,vv,U,L

def fibre_route(vs,areas,target,proof=None):
 am=min(areas)
 if am<=0:return None
 H=E.hull(x for v in vs for x in v)
 if len(H)<3:return None
 def trial(d,N):
  a=hull_setup(H,d)
  if a is None:return False
  n,w,h,vv,U,L=a
  if N==0:return E.PI2*n*am/2*(Q(1,w*w)+Q(1,h*h))>target
  if N not in (16,32,64,128,256):return False
  nodes=sorted(set([Q(i,N) for i in range(N+1)]+[x for x,y in vv])); width=[min(a*x+b for a,b in U)-max(a*x+b for a,b in L) for x in nodes]
  # The polygon's vertical cross-section width is affine BETWEEN consecutive hull x coordinates.
  # Node widths therefore bound every segment exactly; a lower step potential is valid on the entire projection.
  m=Q(target)*2*w*w/(am*n)*(1+Q(1,10**10));phi=(0,0);der=(I,I)
  for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
   step=b-a;W=max(width[i:i+2]);cap=m+4/step**2;v=min(E.PI2/W**2,cap) if W>0 else cap;z=(v-m)*step**2;cs=trans(z)
   if cs is None:return False
   C,S=cs; pn=add(mul(C,phi),mul(mul(enc(step),S),der));dn=add(mul(mul(enc((v-m)*step),S),phi),mul(C,der));phi,der=pn,dn
   if phi[0]<=0:return False
  # All positive nodes, initial derivative 1, oscillatory segment length <=2 < pi (Lemma PS).
  return True
 if proof is not None:
  d=proof['d'];N=proof['N']
  if len(d)!=2 or any(type(x)!=int for x in d) or type(N)!=int:return None
  return proof if trial(tuple(d),N) else None
 # Candidate ordering only. Directions exact integer, selected using rectangle lower bound.
 dd=list(E.dirs(vs));dd.extend(E.vecsub(b,a) for a,b in zip(H,H[1:]+H[:1]));rr=[];seen=set()
 from math import gcd
 for a,b in dd:
  g=gcd(a,b)
  if not g:continue
  d=(a//g,b//g)
  if d[0]<0 or d[0]==0 and d[1]<0:d=(-d[0],-d[1])
  if d in seen:continue
  seen.add(d);n=E.norm2(d);xx=[E.dot(p,d) for p in H];yy=[E.cross(d,p) for p in H];w=max(xx)-min(xx);h=max(yy)-min(yy)
  if not w or not h:continue
  val=E.PI2*n*am/2*(Q(1,w*w)+Q(1,h*h))
  if val>target:return {'d':list(d),'N':0}
  rr.append((val,d))
 for N,top in ((16,4),(32,6),(64,6),(128,6),(256,4)):
  for _,d in sorted(rr,reverse=True)[:top]:
   if trial(d,N):return {'d':list(d),'N':N}
 return None

# Correlated triangle numerator/denominator bound; degree-2 tensor Bernstein products.
P4=list(product((0,1),repeat=4))
@lru_cache(None)
def groups4():
 g={}
 for i,a in enumerate(P4):
  for j,b in enumerate(P4):g.setdefault(tuple(x+y for x,y in zip(a,b)),[]).append((i,j))
 return list(g.values())

def triangle_joint(qs,j,cut,target):
 # All q invariants are multi-affine. A_P*ear_C/A_T proof via degree-elevated products.
 # A_Q = y+D, A_P >= y+D-cut*(D or h) scaled by 2, at fixed upper cut.
 cols=list(zip(*(r[2] for r in qs)));y,h,k,d,m,l=cols;area=[a+b for a,b in zip(y,d)];P=[a-cut*b for a,b in zip(area,d if j==3 else h)]
 for A,B,C in ((h,y,m),(y,k,l),(k,d,tuple(-x for x in m)),(d,h,tuple(-x for x in l))):
  if min(C)<=0:continue
  # TRI*P*C - target*(area*C + A*B)>0 on full box; tensor products exact.
  z=[E.TRI*a-target*b for a,b in zip(P,area)]
  if all(sum(z[i]*C[j]-target*A[i]*B[j] for i,j in g)>0 for g in groups4()):return True
 return False

def label_direct(t,lo,hi,table,proposal=None):
 c,j,lab=t['chart'],t['j'],t['label'];target=E.T
 # First cheap original predicates excluding expensive fibre for O, and original T/TX fibre only if needed later.
 D,s,ll,hh,qs=E.quad(c,list(lo),list(hi));rr=E.ear_ranges(qs)
 if any(b<=0 for a,b in rr[:4]):return {'route':'empty','status':'PROVED'}
 def old():
  r=E.check(c,j,list(lo),list(hi),lab,table)
  return {'route':'old','status':r} if r else None
 # replay routes exactly; proof methods are not accepted merely by their names.
 if proposal is not None:
  route=proposal['route']
  if route=='empty':return None
  if route=='old':return old()
 if lab not in ('T','TX','O'):
  return old() if proposal is None else None
 y,Y=rr[0];h,H=rr[1];d,DM=rr[3];amin=min(r[3] for r in qs)
 cut=Q(hh[4]*hh[5],D*D) if lab=='T' else min(E.untouched(j,rr),default=Q(2))
 tau=min(Q(1),Q(max(0,DM),max(0,y)+max(0,DM))) if j==3 and DM>0 else min(Q(1),Q(H,amin)) if j==0 and amin>0 else Q(1)
 if tau<=0:tau=Q(1)
 den=1-cut*tau
 if lab in ('T','TX'):
  # Strict T/ TX threshold and archived root-table dependency (only std).
  if proposal is None:
   lb=E.lower_quad(c,s,qs,rr)
   if den>0 and lb*den>target:return {'route':'closed','status':'PROVED'}
   if den>0 and c=='std' and table and Q(table)*den>target:return {'route':'table','status':'CONDITIONAL_TABLE'}
  elif proposal['route']=='closed':return {'route':'closed','status':'PROVED'} if den>0 and E.lower_quad(c,s,qs,rr)*den>target else None
  elif proposal['route']=='table':return {'route':'table','status':'CONDITIONAL_TABLE'} if den>0 and c=='std' and table and Q(table)*den>target else None
  if cut<1 and (proposal is None or proposal['route']=='triangle_joint') and triangle_joint(qs,j,cut,target):return {'route':'triangle_joint','status':'PROVED'}
  # q-strip proof of threshold via separate concavity per coordinate, with strict margin.
  if den>0 and (proposal is None or proposal['route']=='q_strip'):
   vs=[((0,0),(s,0),a,b) for a,b,*_ in qs];ars=[a[3] for a in qs]
   if rect_strip(vs,ars,target/den*(1+Q(1,10**10))):return {'route':'q_strip','status':'PROVED'}
  if proposal is not None and proposal['route'] in ('triangle_joint','q_strip'):return None
 # Direct pentagon spectral alternative, valid independently of label; TX may also use X on split child.
 _,vs=E.pent(D,s,ll,hh,qs,j);ars=[E.area2(v) for v in vs]
 if (proposal is None or proposal['route']=='S') and E.strip(vs,ars):return {'route':'S','status':'PROVED'}
 if (proposal is None or proposal['route']=='R') and E.tri_excl(vs,ars):return {'route':'R','status':'PROVED'}
 if lab=='TX' and (proposal is None or proposal['route']=='X') and E.check(c,j,list(lo),list(hi),'X',0):return {'route':'X','status':'PROVED'}
 if proposal is not None and proposal['route'] in ('S','R','X'):return None
 if lab in ('T','TX') and den>0 and (proposal is None or proposal['route']=='q_fibre'):
  qv=[((0,0),(s,0),a,b) for a,b,*_ in qs];qa=[a[3] for a in qs];f=fibre_route(qv,qa,target/den,proposal.get('fibre') if proposal else None)
  if f:return {'route':'q_fibre','status':'PROVED','fibre':f}
 if proposal is None or proposal['route']=='P_fibre':
  f=fibre_route(vs,ars,target,proposal.get('fibre') if proposal else None)
  if f:return {'route':'P_fibre','status':'PROVED','fibre':f}
 return old() if proposal is None else None

def q_direct(t,lo,hi,proposal=None):
 _,s,_,_,qs=E.quad('std',list(lo)+[0,0],list(hi)+[1,1]);rr=E.ear_ranges(qs);L=Q(t['L']);vs=[((0,0),(s,0),a,b) for a,b,*_ in qs];ars=[q[3] for q in qs]
 if any(b<=0 for a,b in rr[:4]):return {'route':'empty','status':'PROVED'}
 if (proposal is None or proposal['route']=='closed') and E.lower_quad('std',s,qs,rr)>=L:return {'route':'closed','status':'PROVED'}
 if (proposal is None or proposal['route']=='q_strip') and rect_strip(vs,ars,L):return {'route':'q_strip','status':'PROVED'}
 if proposal is None or proposal['route']=='q_fibre':
  f=fibre_route(vs,ars,L,proposal.get('fibre') if proposal else None)
  if f:return {'route':'q_fibre','status':'PROVED','fibre':f}
 if (proposal is None or proposal['route']=='legacy_fibre') and E.fibre(vs,ars,L):return {'route':'legacy_fibre','status':'PROVED'}
 return None

def discrete_proof(t,table=0,limit=64):
 direct=(lambda l,h,p=None:q_direct(t,l,h,p)) if t['stage']=='q_other' else (lambda l,h,p=None:label_direct(t,l,h,table,p)); attempts=0
 def rec(lo,hi,b,depth):
  nonlocal attempts
  attempts+=1
  w=direct(lo,hi)
  if w:return {'leaf':w},1,w['status']
  if b<2 or depth>=8 or attempts>=2*limit-1:return None,0,'UNPROVED'
  # Normalize by root widths; alternate dimensions, prioritize geometric w for T/TX if near-equal.
  d=len(lo);axis=max(range(d),key=lambda i: (Q(hi[i])-Q(lo[i]))/(Q(t['hi'][i])-Q(t['lo'][i]))*(Q(5,4) if d==6 and i<4 else 1))
  a,bx=mid_split(lo,hi,axis);l,nl,sl=rec(*a,b-1,depth+1)
  if l is None:return None,0,'UNPROVED'
  r,nr,sr=rec(*bx,b-nl,depth+1)
  if r is None:return None,0,'UNPROVED'
  return {'axis':axis,'left':l,'right':r},nl+nr,'CONDITIONAL_TABLE' if 'CONDITIONAL_TABLE' in (sl,sr) else 'PROVED'
 w,n,s=rec(t['lo'],t['hi'],limit,0)
 dim=len(t['lo'])
 if w is None and limit>=2**dim:
  def uniform(lo,hi,axis):
   nonlocal attempts
   if axis==dim:
    attempts+=1;p=direct(lo,hi)
    return ({'leaf':p},1,p['status']) if p else (None,0,'UNPROVED')
   a,b=mid_split(lo,hi,axis);l,nl,sl=uniform(*a,axis+1)
   if l is None:return None,0,'UNPROVED'
   r,nr,sr=uniform(*b,axis+1)
   if r is None:return None,0,'UNPROVED'
   return {'axis':axis,'left':l,'right':r},nl+nr,'CONDITIONAL_TABLE' if 'CONDITIONAL_TABLE' in (sl,sr) else 'PROVED'
  w,n,s=uniform(t['lo'],t['hi'],0)
 return w,{'parts':n,'attempts':attempts,'status':s}

def discrete_bound(t,w,table=0):
 direct=(lambda l,h,p:q_direct(t,l,h,p)) if t['stage']=='q_other' else (lambda l,h,p:label_direct(t,l,h,table,p))
 def rec(w,lo,hi,depth):
  if depth>8:return 10000,'REJECT'
  if set(w)=={'leaf'}:
   p=direct(lo,hi,w['leaf']);return (1,p['status']) if p==w['leaf'] else (10000,'REJECT')
  if set(w)!={'axis','left','right'} or type(w['axis']) is not int or not 0<=w['axis']<len(lo):return 10000,'REJECT'
  a,b=mid_split(lo,hi,w['axis']);nl,sl=rec(w['left'],*a,depth+1);nr,sr=rec(w['right'],*b,depth+1)
  return nl+nr, 'REJECT' if 'REJECT' in (sl,sr) else 'CONDITIONAL_TABLE' if 'CONDITIONAL_TABLE' in (sl,sr) else 'PROVED'
 n,s=rec(w,t['lo'],t['hi'],0);return s if n<=64 else 'REJECT'
