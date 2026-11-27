#!/usr/bin/env python3
"""Standard-library-only candidate search (slower than optional SciPy engine).
All floating results are proposals; lg_exact is the final proof kernel.
"""
from lg_exact import *
from math import sqrt
from functools import lru_cache
from operator import mul as times

def dot(a,b):return sum(map(times,a,b))
def mv(A,x):return [dot(r,x) for r in A]
def zeros(n,m):return [[0.0]*m for _ in range(n)]
def chol(A):
 n=len(A);L=zeros(n,n)
 for i in range(n):
  for j in range(i+1):
   v=A[i][j]-dot(L[i][:j],L[j][:j]);L[i][j]=sqrt(v) if i==j else v/L[j][j]
 return L
def solve(L,b):
 n=len(b);x=[0.]*n;y=[0.]*n
 for i in range(n):y[i]=(b[i]-dot(L[i][:i],y[:i]))/L[i][i]
 for i in range(n-1,-1,-1):x[i]=(y[i]-sum(L[j][i]*x[j] for j in range(i+1,n)))/L[i][i]
 return x
@lru_cache(None)
def reference(n):
 m=n+1;ms=mons(m);idx={x:i for i,x in enumerate(ms)};H=lcm(*range(2,n+3));M,D=moment(m);size=len(ms)
 def build(N,flux):
  features=[]
  for p in transforms(N,m):
   U=[0]*size
   for i,v in p.items():U[i]=v
   a=[0]*size;b=[0]*size;h=[0]*size;k=[0]*size
   for (i,j),v in zip(ms,U):
    if v:
     if i:a[idx[i-1,j]]=v*i*H
     if j:b[idx[i,j-1]]=v*j*H
     if not flux:h[idx[i+1,j]]=v*(H//(i+j+2));k[idx[i,j+1]]=v*(H//(i+j+2))
   features.append([a,b] if flux else [[x*H for x in U],a,b,h,k])
  return list(zip(*features))
 U=build(n,False);C=build(n+1,True);den=H*H*D
 def gm(a,b):
  tr=[[dot(r,x) for r in M] for x in a];return [[float(Q(dot(x,y),den)) for y in b] for x in tr]
 u,da,db,ah,bh=U;ca,cb=C
 return gm(u,u),gm(da,da),gm(da,db),gm(db,db),gm(ca,ca),gm(ca,cb),gm(cb,cb),gm(ah,ca),gm(bh,cb),gm(bh,ca),gm(ah,cb)
def proposal(chart,lo,hi,j,n,vertices=None,order=1):
 gs,S,_=geometry(chart,lo,hi,j,vertices);k=len(gs);d=len(lo);nu,ui=layout(n,k,True);nc,ci=layout(n+1,k,False);m,aa,ab,bb,ra,rb,rc,xa,xb,xc,xd=reference(n)
 us=[[(i,v) for i,v in enumerate(q) if v>=0] for q in ui];cs=[[(i,v) for i,v in enumerate(q) if v>=0] for q in ci]
 def mat(x):
  K=zeros(nu,nu);M=zeros(nu,nu);R=zeros(nc,nc);B=zeros(nu,nc)
  for t,g in enumerate(gs):
   det,a,b,c=[float(ev(p,x))/S**2 for p in g]
   if det<=0:raise ValueError('nonpositive midpoint determinant')
   for i,ii in us[t]:
    for l,ll in us[t]:K[ii][ll]+=(c*aa[i][l]-b*(ab[i][l]+ab[l][i])+a*bb[i][l])/det;M[ii][ll]+=det*m[i][l]
    for l,ll in cs[t]:B[ii][ll]+=b*(xa[i][l]-xb[i][l])+c*xc[i][l]-a*xd[i][l]
   for i,ii in cs[t]:
    for l,ll in cs[t]:R[ii][ll]+=(c*ra[i][l]-b*(rb[i][l]+rb[l][i])+a*rc[i][l])/det
  return K,M,R,B
 x=[0.]*d;K,M,R,B=mat(x);LK=chol(K);LR=chol(R);u=[1.]*nu
 for _ in range(64):
  v=solve(LK,mv(M,u));c=max(map(abs,v));v=[a/c for a in v];err=max(abs(a-b) for a,b in zip(u,v));u=v
  if err<1e-12:break
 g=mv(M,u);mass=dot(u,g);lam=dot(u,mv(K,u))/mass
 Bt=list(zip(*B));chi=solve(LR,[-x for x in mv(Bt,u)]);uf=[u];cf=[chi]
 if order and d:
  # Shift the first Ritz null direction; orthogonal forcing removes the eigenvalue derivative.
  H=[[K[i][l]-lam*M[i][l]+g[i]*g[l]/mass for l in range(nu)] for i in range(nu)];LH=chol(H)
  def md(A,B):return [[(x-y)/.002 for x,y in zip(a,b)] for a,b in zip(A,B)]
  for i in range(d):
   x[i]=.001;ap=mat(x);x[i]=-.001;am=mat(x);x[i]=0.;Kp,Mp,Rp,Bp=[md(a,b) for a,b in zip(ap,am)];r=[a-lam*b for a,b in zip(mv(Kp,u),mv(Mp,u))];lp=dot(u,r)/mass;up=solve(LH,[lp*a-b for a,b in zip(g,r)]);xp=solve(LR,[-a-b-c for a,b,c in zip(mv(list(zip(*Bp)),u),mv(Bt,up),mv(Rp,chi))]);uf.append(up);cf.append(xp)
 return [[int(round(x*(1<<BITS))) for x in v] for v in uf],[[int(round(x*(1<<BITS))) for x in v] for v in cf]
def choose_mu(p,d):
 a=positive(p,d)
 if a is not None:return [],a
 if not d:return None,None
 sc=max(map(abs,p.values()));g=[p.get(k,0)/sc for k in VARS[:d]];H=[[p.get(VARS[i]+VARS[j],0)/sc*(2 if i==j else 1) for j in range(d)] for i in range(d)]
 def func(mu):
  A=[r.copy() for r in H]
  for i in range(d):A[i][i]+=2*mu[i]
  try:v=solve(chol(A),g);return .5*dot(g,v)+sum(mu),[1-z*z for z in v]
  except (ValueError,ZeroDivisionError):return float('inf'),[0]*d
 mu=[max(0.,(sum(abs(x) for j,x in enumerate(r) if i!=j)-r[i])/2)+1e-8 for i,r in enumerate(H)]
 for _ in range(100):
  v,gr=func(mu);m=[max(0,int(round(x*sc))) for x in mu];b=positive(p,d,m)
  if b is not None:return list(map(str,m)),b
  step=.5
  for i in range(25):
   nxt=[max(0,x-step*y) for x,y in zip(mu,gr)];val,_=func(nxt)
   if val<v:mu=nxt;break
   step/=2
  else:break
 return None,None
def treeproof(ps,d,budget):
 ms=[];failed=[]
 for i,p in enumerate(ps):
  m,b=choose_mu(p,d);ms.append(m)
  if b is None:failed.append(i)
 if not failed:return {'mu':ms},1
 if budget<=1 or not d or any(ps[i].get(0,0)<=0 for i in failed):return None,0
 axis=max(range(d),key=lambda a:sum(abs(v)/max(abs(ps[i].get(0,0)),1) for i in failed for k,v in ps[i].items() if (k>>(5*a))&31))
 a,na=treeproof([split(p,axis,-1) for p in ps],d,budget-1)
 if a is None:return None,0
 b,nb=treeproof([split(p,axis,1) for p in ps],d,budget-na)
 return ({'axis':axis,'left':a,'right':b},na+nb) if b is not None else (None,0)
def certify(chart,lo,hi,j,n,L=TARGET,vertices=None,order=1,budget=64):
 u,chi=proposal(chart,lo,hi,j,n,vertices,order);ds,pd,pc,lb,ub,A,B=assemble(chart,lo,hi,j,n,u,chi,L,vertices,with_base=True);ps=ds+[pd,pc];t,num=treeproof(ps,len(lo),budget)
 if t is None and chart.startswith('q'):
  recorded=Q(L)
  for k in range(1,64):
   L1=recorded*(1+Q(k,16384));p=sub(scale(A,L1.denominator),scale(B,L1.numerator));tmp=ds+[pd,p];t,num=treeproof(tmp,len(lo),1)
   if t is not None:L,pc,ps=L1,p,tmp;break
 info={'lower_center':float(lb) if lb is not None else None,'ritz_center':float(ub) if ub is not None else None,'terms_PD':len(pd),'terms_PC':len(pc),'parts':num,'reason':'positive' if t else 'positivity not established'}
 if t is None:return None,info
 w={'chart':chart,'lo':[str(Q(x)) for x in lo],'hi':[str(Q(x)) for x in hi],'j':j,'N':n,'L':str(L),'u':u,'chi':chi,'digest':digest(ps),'tree':t}
 if vertices is not None:w['vertices']=vertices
 assert replay(w);return w,info
