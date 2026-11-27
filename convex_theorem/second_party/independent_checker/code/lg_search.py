#!/usr/bin/env python3
"""Numerical proposal only. No floating point value is accepted as proof."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import numpy as np
from scipy.linalg import eigh,solve
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
from functools import lru_cache
from math import comb, factorial
from fractions import Fraction as Q
from lg_exact import *
@lru_cache(None)
def ref(n):
 z,w=leggauss(n+3);z=(z+1)/2;w=w/2
 a=(z[:,None]*z[None,:]).ravel();b=(z[:,None]*(1-z[None,:])).ravel();s=a+b;lam=1-s
 weight=(w[:,None]*w[None,:]*z[:,None]).ravel();sw=np.sqrt(weight)[:,None]
 def basis(N):
  return np.stack([factorial(N)/factorial(i)/factorial(j)/factorial(l)*lam**i*a**j*b**l for i,j,l in triples(N)],axis=1)
 def drv(N):
  M=basis(N-1);ind={x:i for i,x in enumerate(triples(N-1))};A=[];B=[];zero=np.zeros(len(a))
  for i,j,l in triples(N):
   t=M[:,ind[i-1,j,l]] if i else zero
   A.append(N*((M[:,ind[i,j-1,l]] if j else zero)-t));B.append(N*((M[:,ind[i,j,l-1]] if l else zero)-t))
  return np.stack(A,axis=1)*sw,np.stack(B,axis=1)*sw
 U=basis(n);idx={x:i for i,x in enumerate(triples(n))};T=np.zeros((len(idx),len(idx)))
 for i,j,l in triples(n):
  k=j+l
  for m in range(k,n+1):
   for J in range(j,j+m-k+1):T[idx[n-m,J,m-J],idx[i,j,l]]=(k+1)/((m+1)*(m+2))*comb(k,j)*comb(m-k,J-j)/comb(m,J)
 h=U@T;da,db=drv(n);ca,cb=drv(n+1);ah=a[:,None]*h*sw;bh=b[:,None]*h*sw;U*=sw
 return U.T@U,da.T@da,da.T@db,db.T@db,ca.T@ca,ca.T@cb,cb.T@cb,ah.T@ca,bh.T@cb,bh.T@ca,ah.T@cb

def proposal(chart,lo,hi,j,n,vertices=None,order=1):
 gs,S,_=geometry(chart,lo,hi,j,vertices); k=len(gs);d=len(lo);nu,ui=layout(n,k,True);nc,ci=layout(n+1,k,False);m,aa,ab,bb,ra,rb,rc,xa,xb,xc,xd=ref(n)
 def geom(x):return [[float(ev(p,x))/S**2 for p in g] for g in gs]
 def mat(x):
  K=np.zeros((nu,nu));M=K.copy();R=np.zeros((nc,nc));B=np.zeros((nu,nc))
  for t,(det,a,b,c) in enumerate(geom(x)):
   if det<=0:raise ValueError('negative midpoint fan determinant')
   us=np.array(ui[t]);cs=np.array(ci[t]);um=us>=0;cm=cs>=0;I=np.ix_(us[um],us[um]);J=np.ix_(cs[cm],cs[cm]);H=np.ix_(us[um],cs[cm]);
   K[I]+=(c*aa-b*(ab+ab.T)+a*bb)[np.ix_(um,um)]/det;M[I]+=det*m[np.ix_(um,um)]
   R[J]+=(c*ra-b*(rb+rb.T)+a*rc)[np.ix_(cm,cm)]/det;B[H]+=(b*(xa-xb)+c*xc-a*xd)[np.ix_(um,cm)]
  return K,M,R,B
 x=[0.]*d;K,M,R,B=mat(x);v,U=eigh(K,M,subset_by_index=[0,0],check_finite=False);u=U[:,0];u=u/np.max(abs(u));chi=solve(R,-B.T@u,assume_a='pos',check_finite=False);us=[u];xs=[chi]
 if order and d:
  g=M@u;H=np.block([[K-v[0]*M,g[:,None]],[g[None,:],np.zeros((1,1))]])
  for i in range(d):
   x[i]=.001;ap=mat(x);x[i]=-.001;am=mat(x);x[i]=0.;Kp,Mp,Rp,Bp=[(a-b)/.002 for a,b in zip(ap,am)]
   up=solve(H,np.r_[-(Kp-v[0]*Mp)@u,0.],assume_a='sym',check_finite=False)[:-1];xp=solve(R,-Bp.T@u-B.T@up-Rp@chi,assume_a='pos',check_finite=False);us.append(up);xs.append(xp)
 return [[int(round(a*(1<<BITS))) for a in row] for row in us],[[int(round(a*(1<<BITS))) for a in row] for row in xs]

def choose_mu(p,d):
 b=positive(p,d)
 if b is not None:return [],b
 # exact kernel first tries absolute bound; empty mu denotes that path
 sc=max((abs(v) for v in p.values()),default=1);g=np.array([p.get(k,0)/sc for k in VARS[:d]]);H=np.array([[p.get(VARS[i]+VARS[j],0)/sc*(2 if i==j else 1) for j in range(d)] for i in range(d)])
 if d==0:return None,None
 tail=sum(abs(v) for k,v in p.items() if sum((k>>(5*i))&31 for i in range(d))>=3)/sc
 if p.get(0,0)/sc<=tail:return None,None
 base=max(0.,-np.linalg.eigvalsh(H)[0]/2)+1e-10
 def fun(mu):
  A=H+np.diag(2*mu)
  try:
   np.linalg.cholesky(A);v=np.linalg.solve(A,g)
   return .5*g@v+sum(mu),1-v*v
  except np.linalg.LinAlgError:return 1e8+sum(abs(mu)),np.ones(d)
 res=minimize(fun,np.full(d,base+max(1e-6,np.max(abs(g))/2)),jac=True,method='L-BFGS-B',bounds=[(0.,None)]*d,options={'maxiter':120,'ftol':1e-13})
 for v in [res.x,np.full(d,base),np.full(d,base+max(1e-6,np.max(abs(g))/2))]:
  mu=[max(0,int(round(float(m)*sc))) for m in v];b=positive(p,d,mu)
  if b is not None:return list(map(str,mu)),b
 return None,None

def treeproof(ps,d,budget):
 def rec(ps,budget):
  ms=[];failed=[]
  for i,p in enumerate(ps):
   m,b=choose_mu(p,d);ms.append(m)
   if b is None:failed.append(i)
  if not failed:return {'mu':ms},1
  if budget<=1 or d==0 or any(ps[i].get(0,0)<=0 for i in failed):return None,0
  scores=[sum(abs(v)/max(abs(ps[i].get(0,0)),1) for k,v in ps[i].items() if (k>>(5*a))&31) for i in failed for a in range(d)]
  axis=max(range(d),key=lambda a:sum(scores[a::d]));a,na=rec([split(p,axis,-1) for p in ps],budget-1)
  if a is None:return None,0
  b,nb=rec([split(p,axis,1) for p in ps],budget-na)
  if b is None:return None,0
  return {'axis':axis,'left':a,'right':b},na+nb
 return rec(ps,budget)

def certify(chart,lo,hi,j,n,L=TARGET,vertices=None,order=1,budget=64):
 u,chi=proposal(chart,lo,hi,j,n,vertices,order);ds,pd,pc,lb,ub,A,B=assemble(chart,lo,hi,j,n,u,chi,L,vertices,with_base=True);ps=ds+[pd,pc]
 t,num=treeproof(ps,len(lo),budget)
 # Nonmonotonic sufficient test: try larger targets if the recorded L is not proved.
 if t is None and chart.startswith('q'):
  recorded=Q(L)
  for k in range(1,64):
   trial=recorded*(1+Q(k,16384));p=sub(scale(A,trial.denominator),scale(B,trial.numerator));family=ds+[pd,p]
   t,num=treeproof(family,len(lo),1)
   if t is not None: L,pc,ps=trial,p,family;break
 info={'lower_center':float(lb) if lb is not None else None,'ritz_center':float(ub) if ub is not None else None,'terms_PD':len(pd),'terms_PC':len(pc),'parts':num,'reason':'positive' if t else 'positivity not established'}
 if t is None:return None,info
 w={'chart':chart,'lo':[str(Q(x)) for x in lo],'hi':[str(Q(x)) for x in hi],'j':j,'N':n,'L':str(L),'u':u,'chi':chi,'digest':digest(ps),'tree':t}
 if vertices is not None:w['vertices']=vertices
 if not replay(w):raise AssertionError('independent kernel rejected candidate')
 return w,info
