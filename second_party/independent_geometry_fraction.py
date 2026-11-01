#!/usr/bin/env python3
"""Original exact-rational interval geometry checker for PAQC mincorner cells.
Only source inputs are frozen JSON numeric literals; NO producer imports.
All operations, spectral root bounds, matrix Loewner tests are Fraction exact.
"""
import math,json,sys,time,itertools
from fractions import Fraction as F
from pathlib import Path
import numpy as np

N=6
SIGNS=list(itertools.product((-1,1),repeat=N))
D2=2**90

def sq_up(x):
    x=F(x)
    if x<0: raise ValueError('negative radicand')
    if x==0:return F(0)
    q=math.isqrt((x.numerator*(D2**2))//x.denominator)
    return F(q+1,D2)
def spectral_max(M):
    a,b,d=M[0][0],(M[0][1]+M[1][0])/2,M[1][1]
    return (a+d)/2+sq_up(((a-d)/2)**2+b*b)
def spectral_min(M):
    a,b,d=M[0][0],(M[0][1]+M[1][0])/2,M[1][1]
    return (a+d)/2-sq_up(((a-d)/2)**2+b*b)
def sym_norm(M):
    a,b,d=M[0][0],(M[0][1]+M[1][0])/2,M[1][1]
    return abs((a+d)/2)+sq_up(((a-d)/2)**2+b*b)
def tr(A):return (A[0][0]+A[1][1])
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def det(A):return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def transpose(A):return [[A[0][0],A[1][0]],[A[0][1],A[1][1]]]
def matmul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def matadd(A,B):return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]
def matsub(A,B):return [[A[i][j]-B[i][j] for j in range(2)] for i in range(2)]
def matscale(A,c):return [[c*A[i][j] for j in range(2)] for i in range(2)]
def inv(A):
 d=det(A)
 if d<=0:raise ValueError('non-PD determinant')
 return [[A[1][1]/d,-A[0][1]/d],[-A[1][0]/d,A[0][0]/d]]
def zero():return [[F(0),F(0)],[F(0),F(0)]]
def genmax(M,G):
 # largest eigenvalue of pencil (M,G) for G>0: exact rational upward formula.
 dg=det(G); dm=det(M)
 b=M[0][0]*G[1][1]+M[1][1]*G[0][0]-2*M[0][1]*G[0][1]
 d=b*b-4*dg*dm
 return (b+sq_up(max(F(0),d)))/(2*dg)

def f64(x):return F(float(x))

def reference():
 theta=2*np.pi*np.arange(5)/5
 XY=np.column_stack((np.cos(theta),np.sin(theta)))
 V=[[f64(t) for t in XY[j]] for j in range(5)]
 C=[sum(v[i] for v in V)/5 for i in range(2)]
 E=[]; a=[]
 for k in range(5):
  p=[V[k][i]-C[i] for i in range(2)];q=[V[(k+1)%5][i]-C[i] for i in range(2)]
  M=[[p[0],q[0]],[p[1],q[1]]]
  E.append(inv(M));a.append(abs(det(M))/2)
 return E,a
E,a0=reference()

def vertices(x):
 return [[F(0),F(0)],[F(1),F(0)]]+[[x[2*i],x[2*i+1]] for i in range(3)]
def Dmats(x):
 V=vertices(x);c=[sum(v[i] for v in V)/5 for i in range(2)];D=[]
 for k in range(5):
  v=[V[k][i]-c[i] for i in range(2)];w=[V[(k+1)%5][i]-c[i] for i in range(2)]
  D.append(matmul([[v[0],w[0]],[v[1],w[1]]],E[k]))
 return D

def Bmats():
 arr=[]
 for i in range(N):
  dx=[F(0)]*6;dx[i]=F(1)
  D1=Dmats(dx); D0=Dmats([F(0)]*6)
  arr.append([matsub(D1[k],D0[k]) for k in range(5)])
 return arr
B=Bmats()

class IV:
 __slots__=('lo','hi')
 def __init__(self,a,b=None):
  if isinstance(a,IV) and b is None:self.lo,self.hi=a.lo,a.hi;return
  self.lo=F(a);self.hi=F(a if b is None else b)
  assert self.lo<=self.hi
 def __add__(self,o):
  o=IV(o);return IV(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return IV(-self.hi,-self.lo)
 def __sub__(self,o):return self+-IV(o)
 def __rsub__(self,o):return IV(o)+-self
 def __mul__(self,o):
  o=IV(o); p=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
  return IV(min(p),max(p))
 __rmul__=__mul__
 def __truediv__(self,o):
  o=IV(o)
  if o.lo<=0<=o.hi:raise ZeroDivisionError
  return self*IV(F(1,o.hi),F(1,o.lo))
 def __rtruediv__(self,o):return IV(o)/self
 def __abs__(self):return IV(0 if self.lo<=0<=self.hi else min(abs(self.lo),abs(self.hi)),max(abs(self.lo),abs(self.hi)))
 def __repr__(self):return f'[{float(self.lo):.10g},{float(self.hi):.10g}]'

def interval_geometry(lo,hi):
 """64 subboxes, directed EXACT rational interval evaluation of 2nd derivatives f=J/A."""
 hw=[hi[i]-lo[i] for i in range(N)]
 all_err=F(0)
 for signs in SIGNS:
  mid=[(lo[i]+hi[i])/2 for i in range(N)]
  x=[IV(lo[i],mid[i]) if signs[i]<0 else IV(mid[i],hi[i]) for i in range(N)]
  Dx=Dmats(x);J=[det(D) for D in Dx]
  if any(z.lo<=0 for z in J):raise ArithmeticError('interval J not positive')
  Aw=sum((IV(a0[k])*J[k] for k in range(5)),IV(0))
  if Aw.lo<=0:raise ArithmeticError('area interval not positive')
  for k in range(5):
   fmin=(J[k]/Aw).lo
   if fmin<=0:raise ArithmeticError('fmin<=0')
   error=F(0)
   for i in range(N):
    Ji=[]
    for z in range(5):
     M=Dx[z];E=B[i][z]
     # derivative of det along rank-one coordinate
     Ji.append(M[1][1]*E[0][0]-M[1][0]*E[0][1]-M[0][1]*E[1][0]+M[0][0]*E[1][1])
    Ai=sum((IV(a0[z])*Ji[z] for z in range(5)),IV(0))
    f2= -IV(2)*Ji[k]*Ai/(Aw*Aw)+IV(2)*J[k]*Ai*Ai/(Aw*Aw*Aw)
    error+=(hw[i]*hw[i]/8)*abs(f2).hi
   all_err=max(all_err,error/fmin)
 return all_err

def check(fname):
 t0=time.time();src=Path('/mnt/data/pent04_workspace/lines/D-topics/probes/PENT/G1')/fname
 info=json.loads(src.read_text());lo=[f64(z) for z in info['cell_lo']];hi=[f64(z) for z in info['cell_hi']]
 Dcorner=[Dmats([hi[i] if sign[i]>0 else lo[i] for i in range(N)]) for sign in SIGNS]
 # convexity all 15 nonadjacent halfplanes, each multiaffine; 64 corners sufficient.
 conv=None
 for sign in SIGNS:
  x=[hi[i] if sign[i]>0 else lo[i] for i in range(N)];V=vertices(x)
  for i in range(5):
   a=V[i];b=V[(i+1)%5]
   for j in range(5):
    if j in (i,(i+1)%5):continue
    c=V[j];o=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    conv=o if conv is None else min(conv,o)
 assert conv>0
 results=[]
 for k in range(5):
  dk=[D[k] for D in Dcorner];J=[det(x) for x in dk]
  assert all(j>0 for j in J) and all(det(B[i][k])==0 for i in range(N))
  X=[inv(x) for x in dk]
  G=[matscale(matmul(xx,transpose(xx)),jj) for xx,jj in zip(X,J)]
  Gcen=matscale([[sum(g[a][b] for g in G) for b in range(2)] for a in range(2)],F(1,64))
  # Find a provable relative Gl=mu Gcen with rational bisection, requiring 64 exact PSD tests.
  LL=F(0);UU=F(1)
  for _ in range(52):
   mi=(LL+UU)/2
   if all(spectral_min(matsub(g,matscale(Gcen,mi)))>=0 for g in G):LL=mi
   else:UU=mi
  mu=LL
  assert mu>0
  Gl=matscale(Gcen,mu)
  assert all(spectral_min(matsub(g,Gl))>=0 for g in G)
  Gi=inv(Gl)
  Rsum=zero();dev=F(0)
  for i in range(N):
   bit=1<<(5-i)
   diffs=[matsub(X[c],X[c^bit]) for c in range(64) if not c&bit]
   R=matscale([[sum(d[a][b] for d in diffs) for b in range(2)] for a in range(2)],F(1,32))
   Rsum=matadd(Rsum,matmul(R,transpose(R)))
   vals=[]
   for dx in diffs:
    e=matsub(dx,R)
    vals.append(spectral_max(matmul(matmul(transpose(e),Gi),e)))
   dev+=max(vals)
  main=genmax(Rsum,Gl)
  best=min((F(1)+tau)*main/F(4)+(F(1)+F(1)/tau)*dev/F(4) for tau in map(F,['0.003','0.01','0.03','0.1','0.3','1.0']))
  eps=max(J)*best
  results.append(dict(eps=eps,mu=mu,gmin=min(spectral_min(g) for g in G)))
  print(fname,'sector',k,'eps',float(eps),'gmin',float(results[-1]['gmin']),'elapsed',round(time.time()-t0,1),flush=True)
 m=1-max(t['eps'] for t in results)
 err=interval_geometry(lo,hi)
 factor=m*(1-err)
 out={'source':fname,'factor_certified_down':math.nextafter(float(factor),-math.inf), 'factor_beats_producer': bool(factor>F(str(info['geometry']['factor']))), 'conditional_LB_beats_Fstar': bool(factor*F.from_float(float(info['Fmin']))>F(str(18.919134195733918))), 'convexity_lower':float(conv),'sectors':[{k:float(v) for k,v in x.items()} for x in results], 'm':float(m),'interp_err':float(err),'factor_lower':float(factor),'producer_factor':info['geometry']['factor'],'producer_Fmin':info['Fmin'],'derived_LB':float(factor*F(str(info['Fmin']))),'producer_LB':info['LB'],'Fstar_up':18.919134195733918,'sec':time.time()-t0}
 print('RESULT',json.dumps(out,indent=2),flush=True)
 return out
if __name__=='__main__':
 output=check(sys.argv[1]);Path('/mnt/data/PENT-04-GPT-evidence/logs/geom_'+sys.argv[1]).write_text(json.dumps(output,indent=2))

