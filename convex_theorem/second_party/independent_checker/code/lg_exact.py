#!/usr/bin/env python3
"""Independent LG proof kernel. Standard library only. Floating-point search is separate."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial as f, comb, lcm
from itertools import product
import hashlib,json,operator
KS=Q('36.33682907'); TARGET=Q('18.919134195733914'); BITS=32
VARS=tuple(1<<(5*i) for i in range(6))
def add(*ps):
 d={}
 for p in ps:
  for k,v in p.items(): d[k]=d.get(k,0)+v
 return {k:v for k,v in d.items() if v}
def scale(p,c): return {k:v*c for k,v in p.items() if v*c}
def mul(p,q):
 if len(p)>len(q): p,q=q,p
 d={}
 for k,v in p.items():
  for j,w in q.items(): d[k+j]=d.get(k+j,0)+v*w
 return {k:v for k,v in d.items() if v}
def multi(ps):
 d={0:1}
 for p in ps: d=mul(d,p)
 return d
def sub(p,q): return add(p,scale(q,-1))
def ev(p,x): return sum(v*prod(x[i]**((k>>(5*i))&31) for i in range(len(x))) for k,v in p.items())
def prod(xs):
 a=1
 for x in xs:a*=x
 return a
@lru_cache(None)
def triples(n): return tuple((i,j,n-i-j) for i in range(n+1) for j in range(n-i+1))
@lru_cache(None)
def layout(n,k,bdy):
 ids={}; rows=[]
 for t in range(k):
  row=[]
  for i,j,l in triples(n):
   if bdy and i==0: row.append(-1);continue
   if i==n:
    if not bdy: row.append(-1);continue # gauge chi(O)=0
    key=('center',)
   elif j==0: key=('spoke',(t+1)%k,i)
   elif l==0: key=('spoke',t,i)
   elif i==0:key=('edge',t,j)
   else:key=('interior',t,i,j)
   if key not in ids:ids[key]=len(ids)
   row.append(ids[key])
  rows.append(row)
 return len(ids),rows
@lru_cache(None)
def mons(n):return tuple((a,d-a) for d in range(n+1) for a in range(d+1))
@lru_cache(None)
def transforms(n,m):
 idx={x:i for i,x in enumerate(mons(m))}; rows=[]
 for i,j,l in triples(n):
  d={}; c=f(n)//(f(i)*f(j)*f(l))
  for p in range(i+1):
   for q in range(i-p+1): d[idx[j+p,l+q]]=c*(f(i)//(f(p)*f(q)*f(i-p-q)))*(-1)**(p+q)
  rows.append(d)
 return rows
@lru_cache(None)
def moment(m):
 ms=mons(m); D=f(2*m+2)
 return tuple(tuple(f(a+c)*f(b+d)*(D//f(a+b+c+d+2)) for c,d in ms) for a,b in ms),D

def geometry(chart,lo,hi,j=3,vertices=None,origin=None):
 d=len(lo); x=[{0:(Q(a)+Q(b))/2,VARS[i]:(Q(b)-Q(a))/2} for i,(a,b) in enumerate(zip(lo,hi))]
 z={};one={0:Q(1)}
 if vertices is not None: vs=[({0:Q(a)},{0:Q(b)}) for a,b in vertices]
 elif chart=='near':vs=[(z,z),(one,z),*[(x[i],x[i+1]) for i in (0,2,4)]]
 else:
  if chart in ('std','qstd'):q2=(x[0],x[1]);q3=(x[2],x[3])
  elif chart in ('E2','qE2'):q3=(x[0],x[1]);q2=(add(one,mul(x[2],sub(q3[0],one))),add(mul(x[2],q3[1]),x[3]))
  elif chart in ('E3','qE3'):q2=(x[0],x[1]);q3=(mul(x[2],q2[0]),add(mul(x[2],q2[1]),x[3]))
  else:raise ValueError(chart)
  vs=[(z,z),(one,z),q2,q3]
  if not chart.startswith('q'):
   t,u=x[4:6]
   if j==3:vs=[*vs[:3],tuple(add(a,mul(t,sub(b,a))) for a,b in zip(q3,q2)),tuple(mul(sub(one,u),a) for a in q3)]
   elif j==0:vs=[(t,z),*vs[1:],tuple(mul(u,a) for a in q3)]
   else:raise ValueError(j)
 o=[Q(sum(v[i].get(0,0) for v in vs))/len(vs) for i in range(2)] if origin is None else list(map(Q,origin))
 e=[tuple(sub(p,{0:o[i]}) for i,p in enumerate(v)) for v in vs]
 S=lcm(*(Q(v).denominator for q in e for p in q for v in p.values()))
 e=[tuple({k:int(v*S) for k,v in p.items()} for p in q) for q in e]
 gs=[]
 for a,b in zip(e,e[1:]+e[:1]):
  dot=lambda a,b:add(mul(a[0],b[0]),mul(a[1],b[1]))
  gs.append((sub(mul(a[0],b[1]),mul(a[1],b[0])),dot(a,a),dot(a,b),dot(b,b)))
 return gs,S,vs

def feature_rows(n,k,u,chi):
 m=n+1; ms=mons(m); idx={x:i for i,x in enumerate(ms)}; H=lcm(*range(2,n+3)); nu,ui=layout(n,k,True); nc,ci=layout(n+1,k,False)
 if len(u)!=len(chi) or not 1<=len(u)<=7 or any(len(x)!=nu for x in u) or any(len(x)!=nc for x in chi):raise ValueError('coefficient dimensions')
 if not all(type(x)==int for row in u+chi for x in row):raise ValueError('coefficients not integer')
 for t in range(k):
  rows=[]
  for uu,cc in zip(u,chi):
   def mono(vals,ids,N):
    a=[0]*len(ms)
    for ix,p in zip(ids,transforms(N,m)):
     if ix>=0 and vals[ix]:
      for j,w in p.items():a[j]+=w*vals[ix]
    return a
   U=mono(uu,ui[t],n); C=mono(cc,ci[t],n+1)
   def deriv(a,d):
    out=[0]*len(ms)
    for (p,q),v in zip(ms,a):
     if d==0 and p:out[idx[p-1,q]]=v*p*H
     if d==1 and q:out[idx[p,q-1]]=v*q*H
    return out
   def hm(d):
    out=[0]*len(ms)
    for (p,q),v in zip(ms,U):
     if v:out[idx[p+(d==0),q+(d==1)]]=v*(H//(p+q+2))
    return out
   rows.append([[x*H for x in U],deriv(U,0),deriv(U,1),hm(0),hm(1),deriv(C,0),deriv(C,1)])
  yield rows

def forms(n,k,u,chi):
 M,D=moment(n+1); keys=(0,)+VARS[:len(u)-1]; dot=lambda a,b:sum(map(operator.mul,a,b))
 for rows in feature_rows(n,k,u,chi):
  # row features are exact integer coefficients over common denominator H*2^BITS.
  tr=[[[dot(a,v) for a in M] for v in r] for r in rows]
  out=[]
  for a,b in [(0,0),(1,1),(1,2),(2,2),(3,3),(3,4),(4,4),(3,5),(4,6),(4,5),(3,6),(5,5),(5,6),(6,6)]:
   p={}
   for i,ki in enumerate(keys):
    for j,kj in enumerate(keys):
     key=ki+kj;p[key]=p.get(key,0)+dot(tr[i][a],rows[j][b])
   out.append({key:v for key,v in p.items() if v})
  yield out

def assemble(chart,lo,hi,j,n,u,chi,L,vertices=None,origin=None,with_base=False):
 gs,S,_=geometry(chart,lo,hi,j,vertices,origin);k=len(gs); ds=[g[0] for g in gs]; dp=multi(ds); sd=add(*ds)
 A1=[]; A0=[]; A2q=[]; A2x=[]; A2r=[]
 for i,((det,a,b,c),fs) in enumerate(zip(gs,forms(n,k,u,chi))):
  mass,aa,ab,bb,qaa,qab,qbb,xa1,xb2,xb1,xa2,raa,rab,rbb=fs
  others=multi(ds[:i]+ds[i+1:])
  metric=lambda aa,ab,bb:add(mul(c,aa),scale(mul(b,ab),-2),mul(a,bb))
  A1.append(mul(det,mass)); A0.append(mul(metric(aa,ab,bb),others))
  A2q.append(mul(det,add(mul(a,qaa),scale(mul(b,qab),2),mul(c,qbb))))
  A2x.append(add(mul(b,sub(xa1,xb2)),mul(c,xb1),scale(mul(a,xa2),-1)))
  A2r.append(mul(metric(raa,rab,rbb),others))
 a1=add(*A1);a0=add(*A0);a2=add(mul(dp,add(add(*A2q),scale(add(*A2x),2*S*S))),scale(add(*A2r),S**4))
 K,J=KS.numerator,KS.denominator; L=Q(L); F,G=L.numerator,L.denominator
 a1d=mul(a1,dp);pd=sub(scale(a1d,2*K),scale(mul(sd,a0),J));baseA=mul(sd,pd);baseB=scale(sub(scale(a2,2*K),scale(mul(sd,a1d),J)),2);pc=sub(scale(baseA,G),scale(baseB,F))
 # Independent exact center evaluation from reference integrals.
 area=Q(sd.get(0,0),2*S*S); H=lcm(*range(2,n+3));Den=(H*(1<<BITS))**2*moment(n+1)[1]
 d0=dp.get(0,0)
 if d0:
  b1=Q(a1.get(0,0),Den*S*S); b0=Q(a0.get(0,0),Den*d0); b2=Q(a2.get(0,0),Den*S**4*d0)
  dc=KS*b1-area*b0; dc2=KS*b2-area*b1
  lower=area*dc/dc2 if dc>0 and dc2>0 else None
  upper=area*b0/b1 if b1 else None
 else:lower=upper=None
 return (ds,pd,pc,lower,upper,baseA,baseB) if with_base else (ds,pd,pc,lower,upper)

def positive(p,d,mu=None):
 """Return exact lower bound, or None. Quadratic dual mu is rational, nonnegative."""
 c=p.get(0,0); a=c-sum(abs(v) for k,v in p.items() if k)
 if a>0:return Q(a)
 if mu is None:return None
 if len(mu)!=d or any(Q(x)<0 for x in mu):return None
 g=[Q(p.get(k,0)) for k in VARS[:d]]; H=[[Q(p.get(VARS[i]+VARS[j],0))*(2 if i==j else 1) for j in range(d)] for i in range(d)]
 tail=sum(abs(v) for k,v in p.items() if sum((k>>(5*i))&31 for i in range(d))>=3)
 for i in range(d):H[i][i]+=2*Q(mu[i])
 # Exact elimination computes g^T H^-1 g; reject nonpositive pivots.
 s=Q(0)
 for i in range(d):
  a=H[i][i]
  if a<=0:return None
  s+=g[i]**2/a
  for j in range(i+1,d):
   t=H[j][i]/a;g[j]-=t*g[i]
   for k in range(i+1,d):H[j][k]-=t*H[i][k]
 lb=c-tail-s/2-sum(map(Q,mu))
 return lb if lb>0 else None

def split(p,axis,sign):
 shift=5*axis; unit=1<<shift; M=max(((k>>shift)&31 for k in p),default=0);r={}
 for k,v in p.items():
  a=(k>>shift)&31; k-=a*unit;v<<=M-a
  for j in range(a+1):
   key=k+j*unit;r[key]=r.get(key,0)+v*comb(a,j)*sign**(a-j)
 return {k:v for k,v in r.items() if v}
def digest(ps):
 h=hashlib.sha256()
 for p in ps:
  for k,v in sorted(p.items()):h.update(f'{k}:{v}\n'.encode())
  h.update(b'END\n')
 return h.hexdigest()
def replay(w):
 lo,hi=w['lo'],w['hi']; d=len(lo)
 if d>6 or len(hi)!=d or any(Q(a)>=Q(b) for a,b in zip(lo,hi)):return False
 if len(w['u']) not in (1,d+1):return False
 ds,pd,pc,*_=assemble(w['chart'],lo,hi,w['j'],w['N'],w['u'],w['chi'],w['L'],w.get('vertices'))
 ps=ds+[pd,pc]
 if digest(ps)!=w['digest']:return False
 # det-positive at every point; geometry winding is separately established by application scope.
 def walk(ps,t):
  if 'axis' in t:
   i=t['axis'];return (walk([split(p,i,-1) for p in ps],t['left'])+walk([split(p,i,1) for p in ps],t['right'])) if 0<=i<d else 1000
  if len(t['mu'])!=len(ps) or any(positive(p,d,m) is None for p,m in zip(ps,t['mu'])):return 1000
  return 1
 return walk(ps,w['tree'])<=64
