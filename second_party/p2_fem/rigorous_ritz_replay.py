#!/usr/bin/env python3
"""Independent strict conforming P2 variational upper bound on regular pentagon.
Uses dyadic-rational mesh/trial data + directed-rounded Decimal interval arithmetic.
Only reads independently generated Ritz trial .npz, NEVER the bundled code.
Proof: Q (rational pentagon) strictly inside exact regular P5; trial belongs to H0^1(Q),
λ1(P5) ≤ λ1(Q) ≤ a(w,w)/b(w,w); multiply A(P5) by interval bound.
"""
import numpy as np
from fractions import Fraction as F
from decimal import Decimal as D,Context,ROUND_FLOOR,ROUND_CEILING,localcontext
from pathlib import Path
import math,sys,json,time
from mpmath import iv
iv.dps=55
LF=Context(prec=70,rounding=ROUND_FLOOR); UF=Context(prec=70,rounding=ROUND_CEILING)

class I:
 __slots__=('lo','hi')
 def __init__(self,lo,hi=None):self.lo=D(lo);self.hi=self.lo if hi is None else D(hi)
 @staticmethod
 def fr(z):
  if not isinstance(z,F):z=F(z)
  return I(LF.divide(D(z.numerator),D(z.denominator)),UF.divide(D(z.numerator),D(z.denominator)))
 def __add__(self,o):
  if not isinstance(o,I):o=I.fr(o)
  return I(LF.add(self.lo,o.lo),UF.add(self.hi,o.hi))
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+(-o if isinstance(o,I) else -I.fr(o))
 def __rsub__(self,o):return I.fr(o)+(-self)
 def __mul__(self,o):
  if not isinstance(o,I):o=I.fr(o)
  vals=[(x,y) for x in (self.lo,self.hi) for y in (o.lo,o.hi)]
  return I(min(LF.multiply(x,y) for x,y in vals),max(UF.multiply(x,y) for x,y in vals))
 __rmul__=__mul__
 def __truediv__(self,o):
  if not isinstance(o,I):o=I.fr(o)
  assert o.lo>0 or o.hi<0
  return self*I(min(LF.divide(D(1),o.lo),LF.divide(D(1),o.hi)),max(UF.divide(D(1),o.lo),UF.divide(D(1),o.hi)))
 def __repr__(self):return f'[{self.lo}, {self.hi}]'

# Algebraically exact simplex reference moments for monomials (x^i y^j).
EXP=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
def integ(i,j):return F(math.factorial(i)*math.factorial(j),math.factorial(i+j+2))
MOM=[[I.fr(integ(EXP[a][0]+EXP[b][0],EXP[a][1]+EXP[b][1])) for b in range(6)] for a in range(6)]
T3=[[I.fr(integ(a,b)) for a,b in [(0,0),(1,0),(0,1)]],]

def integral_square(coeff):
 res=I(0)
 for a in range(len(coeff)):
  for b in range(a,len(coeff)):
   i=EXP[a][0]+EXP[b][0] if len(coeff)==6 else ([0,1,0][a]+[0,1,0][b])
   j=EXP[a][1]+EXP[b][1] if len(coeff)==6 else ([0,0,1][a]+[0,0,1][b])
   res=res+I.fr((2 if a!=b else 1)*integ(i,j))*coeff[a]*coeff[b]
 return res

def calc_trial(npz,m):
 a=np.load(npz);V=a['V'];fv=[(F.from_float(float(p[0])),F.from_float(float(p[1]))) for p in V]
 # reconstruct topology and genuinely exact points from barycentric rational mesh
 idx={};coords={};tri=[];nverts=0
 def put(k,aa,bb):
  p=(aa*V[k]+bb*V[(k+1)%5])/m;key=tuple(np.round(p,14))
  px=((aa*fv[k][0]+bb*fv[(k+1)%5][0])/m,(aa*fv[k][1]+bb*fv[(k+1)%5][1])/m)
  if key not in idx:idx[key]=len(idx);coords[idx[key]]=px
  else: assert coords[idx[key]]==px
  return idx[key]
 for k in range(5):
  t={}
  for aa in range(m+1):
   for bb in range(m+1-aa):t[(aa,bb)]=put(k,aa,bb)
  for aa in range(m):
   for bb in range(m-aa):
    tri.append([t[(aa,bb)],t[(aa+1,bb)],t[(aa,bb+1)]])
    if aa+bb<m-1:tri.append([t[(aa+1,bb)],t[(aa+1,bb+1)],t[(aa,bb+1)]])
 assert len(idx)==len(a['pts'])-len(set((min(int(t[i]),int(t[j])),max(int(t[i]),int(t[j]))) for t in tri for i,j in ((0,1),(1,2),(2,0))))
 assert len(tri)==len(a['els'])
 edges={}
 for tri0 in tri:
  mids=[]
  for i,j in ((0,1),(1,2),(2,0)):
   ea,eb=sorted((tri0[i],tri0[j]));key=(ea,eb)
   if key not in edges:
    mid=len(idx)+len(edges);edges[key]=mid
    coords[mid]=((coords[ea][0]+coords[eb][0])/2,(coords[ea][1]+coords[eb][1])/2)
   mids.append(edges[key])
  # exact triangle vertex ordering MUST equal stored connectivity
  if 'elements' not in locals():elements=[]
  elements.append([*tri0,*mids])
 assert np.array_equal(np.array(elements),a['els'])
 assert len(coords)==len(a['pts'])
 for i,p in coords.items():
  assert abs(float(p[0])-a['pts'][i,0])<1e-12 and abs(float(p[1])-a['pts'][i,1])<1e-12
 vtrial=[F(0)]*len(coords)
 for i,t in zip(a['ids'],a['trial']):vtrial[int(i)]=F.from_float(float(t))
 assert all(vtrial[int(i)]==0 for i in np.where(a['boundary'])[0])
 # exact zero trace: boundary mesh nodes lie on polygon sides and their coefficients = 0
 NUM=I(0);DEN=I(0)
 for n,t in enumerate(elements):
  A,B,C=[coords[i] for i in t[:3]]
  bx,by=B[0]-A[0],B[1]-A[1];cx,cy=C[0]-A[0],C[1]-A[1]
  det=bx*cy-by*cx
  assert det>0
  ts=[I.fr(vtrial[i]) for i in t]
  c0=ts[0]
  c1=-3*ts[0]-ts[1]+4*ts[3]
  c2=-3*ts[0]-ts[2]+4*ts[5]
  c3=2*ts[0]+2*ts[1]-4*ts[3]
  c4=4*ts[0]-4*ts[3]+4*ts[4]-4*ts[5]
  c5=2*ts[0]+2*ts[2]-4*ts[5]
  coeff=[c0,c1,c2,c3,c4,c5]
  # dxi = (c1, 2c3, c4); deta=(c2,c4,2c5)
  dxi=[c1,2*c3,c4];deta=[c2,c4,2*c5]
  g1=[I.fr(cy)*dxi[i]-I.fr(by)*deta[i] for i in range(3)]
  g2=[I.fr(-cx)*dxi[i]+I.fr(bx)*deta[i] for i in range(3)]
  NUM=NUM+(integral_square(g1)+integral_square(g2))/I.fr(det)
  DEN=DEN+I.fr(det)*integral_square(coeff)
  if n and n%2500==0:print('completed elements',n,'/',len(elements),flush=True)
 assert NUM.lo>0 and DEN.lo>0
 lam=NUM/DEN
 return fv,lam,NUM,DEN,len(elements)

def ivfrac(z):return iv.mpf(z.numerator)/iv.mpf(z.denominator)
def verify_inclusion(fv):
 pi=iv.pi; W=[(iv.cos(2*pi*i/5),iv.sin(2*pi*i/5)) for i in range(5)]
 positive=[]
 for i in range(5):
  x,y=W[i];xx,yy=W[(i+1)%5]
  for px,py in fv:
   cross=(xx-x)*(ivfrac(py)-y)-(yy-y)*(ivfrac(px)-x)
   positive.append(float(cross.a))
   assert cross.a>0,('Q not inside target',i,px,py,str(cross))
 area=iv.mpf(5)/2*iv.sin(2*pi/5)
 return area,min(positive)

if __name__=='__main__':
 npz=sys.argv[1];m=int(sys.argv[2]);out=sys.argv[3]
 t0=time.time();fv,lam,num,den,ne=calc_trial(npz,m);area,margin=verify_inclusion(fv)
 # iv arithmetic on lambda interval from Decimal interval endpoints; these are decimal exact
 lam_iv=iv.mpf([str(lam.lo),str(lam.hi)]);Flam=lam_iv*area
 res=dict(m=m,triangles=ne,area_reg_interval=str(area),P5_contains_Q=True,min_side_cross=margin,
          numerator_interval=repr(num),denominator_interval=repr(den),lambda_upper=str(lam.hi),F_interval=str(Flam),
          Fstar_up_required='18.919134195733918',PASS=Flam.b < iv.mpf('18.919134195733918'),seconds=round(time.time()-t0,2))
 Path(out).write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))

