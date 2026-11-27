#!/usr/bin/env python3
"""Independent far label sufficient tests from FAR_PROOFS. Exact integer/Fraction comparisons.
T/TX may be conditional on archived table LB. O uses rectangle then rigorous fibre transfer.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt,cos,sin,pi,lcm
from lg_exact import TARGET
# Machin identity, alternating partial sums with even number of terms, lower/upper.
def atan_bounds(k,n=24):
 s=sum((Q((-1)**i,(2*i+1)*k**(2*i+1)) for i in range(n)),Q(0));t=s+Q((-1)**n,(2*n+1)*k**(2*n+1));return min(s,t),max(s,t)
a,b=atan_bounds(5);c,d=atan_bounds(239);PI_LOW=16*a-4*d; PI2=(PI_LOW**2); PI2=Q((PI2*10**30).__floor__(),10**30)
TRI=4*PI2/Q(1732050807568877293527446341506,10**30) # upper sqrt3
assert Q(1732050807568877293527446341506,10**30)**2>3
T=TARGET
C=[Q(x) for x in [1.3090169943749475,.9510565162951536,.5,1.5388417685876268,1-1.3090169943749475,.9510565162951536]]
def vecsub(a,b):return (a[0]-b[0],a[1]-b[1])
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def norm2(a):return dot(a,a)
def rng(a):
 a=tuple(a);return min(a),max(a)
def minabs(a,b):return 0 if a<=0<=b else min(abs(a),abs(b))
def sqrtup(f,bits=40):
 f=Q(f);m=1<<bits;n=isqrt(f.numerator*m*m//f.denominator);return Q(n+(n*n*f.denominator<f.numerator*m*m),m)
def sqrtlo(f,bits=40):
 f=Q(f);m=1<<bits;return Q(isqrt(f.numerator*m*m//f.denominator),m)
def quad(chart,lo,hi):
 qq=list(map(Q,lo+hi));D=lcm(*(x.denominator for x in qq));s=D*D;ll=[int(x*D) for x in qq[:6]];hh=[int(x*D) for x in qq[6:]];qs=[]
 for a,b,c,d in product(*zip(ll[:4],hh[:4])):
  if chart=='std':q2=(a*D,b*D);q3=(c*D,d*D)
  elif chart=='E2':q3=(a*D,b*D);q2=(s+c*(a-D),c*b+d*D)
  elif chart=='E3':q2=(a*D,b*D);q3=(c*a,c*b+d*D)
  else:raise ValueError(chart)
  x,y=q2;z,h=q3;dd=cross(q2,q3);k=dd+(y-h)*s;m=dd-h*s;l=(h-y)*s
  qs.append((q2,q3,(y*s,h*s,k,dd,m,l),y*s+dd))
 return D,s,ll,hh,qs

def pent(D,s,ll,hh,qs,j):
 V=[]; S=s*D
 for q2,q3,_,_ in qs:
  a,b=q2;c,d=q3;A=(a*D,b*D);B=(c*D,d*D)
  for t,u in product(*zip(ll[4:],hh[4:])):
   if j==3:V.append(((0,0),(S,0),A,(c*D+t*(a-c),d*D+t*(b-d)),((D-u)*c,(D-u)*d)))
   else:V.append(((t*s,0),(S,0),A,B,(u*c,u*d)))
 return S,V

def ear_ranges(qs):return [rng(v) for v in zip(*(q[2] for q in qs))]
def untouched(j,rr):
 (y,Y),(h,H),(k,K),(d,D),(m,M),(l,L)=rr;y=max(0,y);h=max(0,h);k=max(0,k);d=max(0,d)
 pairs=[(Y*H,d*m),(Y*K,d*l)] if j==3 else [(Y*K,h*l),(D*K,h*(-M))]
 cond=[d>0 and m>0,d>0 and l>0] if j==3 else [h>0 and l>0,h>0 and M<0]
 return [Q(max(0,a),b) for (a,b),ok in zip(pairs,cond) if ok and b>0]
def stub(j,D,ll,hh,qs):
 ok=[True,True]
 for _,_,(y,h,k,d,m,l),_ in qs:
  for t,u in product(*zip(ll[4:],hh[4:])):
   if j==3:
    C1=t*k-u*m;P1=k*(D-2*t)*D+t*u*m;C2=u*h-t*l;P2=h*(D-2*u)*D+t*u*l
   else:
    C1=u*m+t*y;P1=y*(D-2*t)*D-t*u*m;C2=u*d-t*l;P2=d*(D-2*u)*D+t*u*l
   ok[0]&=C1>0 and P1<0;ok[1]&=C2>0 and P2<0
   if not any(ok):return False
 return any(ok)
# interval rational operations for E3 homogeneous face. r<=1 and s=1/r<=1.
def ia(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def im(a,b):return rng([x*y for x in a for y in b])
def iv(c):return (Q(c),Q(c))
def homface(lo,hi,j,tu):
 if lo[1]!=0 or lo[3]!=0 or lo[2]<0:return False
 x=tuple(map(Q,[lo[0],hi[0]]));b=tuple(map(Q,[lo[2],hi[2]]));one=iv(1);xm=ia(x,neg(one));bm=ia(one,neg(b))
 def onepiece(r,reverse):
  if reverse:y=r;h=ia(im(b,r),one);k=ia(xm,im(bm,r));d=x;l=ia(one,neg(im(bm,r)));m=ia(xm,neg(im(b,r)))
  else:y=one;h=ia(b,r);k=ia(im(r,xm),bm);d=im(x,r);l=ia(r,neg(bm));m=ia(im(r,xm),neg(b))
  if k[1]<=0:return True
  pairs=[(im(y,h),im(d,m),m),(im(y,k),im(d,l),l)] if j==3 else [(im(y,k),im(h,l),l),(im(d,k),im(h,neg(m)),neg(m))]
  return any(sig[0]>0 and tu*den[0]>num[1] for num,den,sig in pairs)
 # Endpoint degeneracies may defeat this conservative interval test; report unproved.
 return all(onepiece((Q(i,64),Q(i+1,64)),rev) for rev in (False,True) for i in range(64))
def lower_quad(chart,s,qs,rr):
 y,Y=rr[0];h,H=rr[1];k,K=rr[2];d,D=rr[3];m,M=rr[4];l,L=rr[5]; amin,amax=rng(q[3] for q in qs);w=max(Y,H);lb=2*PI2
 if w>0:lb=max(lb,PI2*s*s/(2*w),PI2*max(0,amin)*s*s/(2*w*w))
 if amin>0:
  for aa,bb,cc in [(H,Y,m),(Y,K,l),(K,D,-M),(D,H,-L)]:
   if cc>0:lb=max(lb,TRI/(1+Q(max(0,aa)*max(0,bb),cc*amin)))
 if chart=='std' and amax>0:
  edges=[[(q[0][0]-s,q[0][1]) for q in qs],[(q[1][0]-q[0][0],q[1][1]-q[0][1]) for q in qs],[q[1] for q in qs]]
  per=1+sum(sqrtlo(Q(sum(minabs(*rng(col))**2 for col in zip(*e)),s*s)) for e in edges);lb=max(lb,PI2*per*per*s*s/(8*amax))
 return lb

def area2(v):return sum(cross(a,b) for a,b in zip(v,v[1:]+v[:1]))
def geom_N(vs):
 # exact mean of corner vertices. direction is parallel to midpoint edge.
 mid=[(sum(v[i][0] for v in vs),sum(v[i][1] for v in vs)) for i in range(5)]
 for rev in (1,-1):
  for start in range(5):
   ids=[(start+rev*i)%5 for i in range(5)];a,b=ids[:2];e=vecsub(mid[b],mid[a]);ee=norm2(e)
   if ee==0:continue
   ds=[vecsub(v[b],v[a]) for v in vs];md=min(dot(z,e) for z in ds)
   if md<=0:continue
   total=Q(0)
   for p,ii in enumerate(ids[2:]):
    z=vecsub(mid[ii],mid[a]);rx=dot(z,e);ry=cross(e,z)
    sm=0
    for v,(dx,dy) in zip(vs,ds):
     nx,ny=vecsub(v[ii],v[a]);sm=max(sm,(ee*nx-rx*dx+ry*dy)**2+(ee*ny-rx*dy-ry*dx)**2)
    rho=Q(isqrt(sm*ee)+(isqrt(sm*ee)**2<sm*ee),ee*md)
    delt=sqrtup((Q(rx,ee)-C[2*p])**2+(Q(rev*ry,ee)-C[2*p+1])**2);total+=(rho+delt+Q(1,1<<48))**2
    if total>=Q(9,25):break
   else:return True
 return False
FIXED=tuple((round(1024*cos(i*pi/36)),round(1024*sin(i*pi/36))) for i in range(36))
def dirs(vs):
 mid=[(sum(v[i][0] for v in vs),sum(v[i][1] for v in vs)) for i in range(len(vs[0]))]
 return [(-d[1],d[0]) for a,b in __import__('itertools').combinations(mid,2) if (d:=vecsub(a,b))!=(0,0)]+list(FIXED)
def strip(vs,areas):
 # edge strip first; no positivity of signed corner areas assumed.
 p,q=PI2.numerator*T.denominator,PI2.denominator*T.numerator
 for i in range(5):
  es=[vecsub(v[(i+1)%5],v[i]) for v in vs];l2=sum(minabs(*rng(col))**2 for col in zip(*es));c=max(cross(e,vecsub(w,v[i])) for v,e in zip(vs,es) for w in v)
  if c>0 and p*l2>2*q*c:return True
 if min(areas)<=0:return False
 for n in dirs(vs):
  d2=norm2(n)
  if d2 and all(p*d2*a>2*q*(max(pr)-min(pr))**2 for a,v in zip(areas,vs) for pr in [[dot(n,x) for x in v]]):return True
 return False

def tri_excl(vs,areas):
 am=min(areas)
 if am<=0:return False
 es=[[vecsub(v[(i+1)%5],v[i]) for i in range(5)] for v in vs];ears=[]
 for i in range(5):
  A=max(cross(e[(i-1)%5],e[i]) for e in es);B=max(cross(e[i],e[(i+1)%5]) for e in es);C0=min(cross(e[(i-1)%5],e[(i+1)%5]) for e in es)
  ears.append(Q(A*B,C0) if C0>0 and A>0 and B>0 else None)
 return any(a is not None and b is not None and TRI*am>T*(am+a+b) for i,a in enumerate(ears) for b in [ears[(i+2)%5]])

def hull(ps):
 ps=sorted(set(ps));a=[];b=[]
 for seq,out in [(ps,a),(list(reversed(ps)),b)]:
  for x in seq:
   while len(out)>1 and cross(vecsub(out[-1],out[-2]),vecsub(x,out[-1]))<=0:out.pop()
   out.append(x)
 return a[:-1]+b[:-1]
# Directed fixed-point interval arithmetic for transfer in scaled independent variable.
I=1<<48
def enc(x):
 x=Q(x)*I;return x.__floor__(),x.__ceil__()
def iadd(a,b):return a[0]+b[0],a[1]+b[1]
def imul(a,b):
 lo,hi=rng(x*y for x in a for y in b);return lo//I,-((-hi)//I)
def idiv(a,k):return a[0]//k,-((-a[1])//k)
def transfer(z):
 # C(z)=sum z^k/(2k)!, S(z)=sum z^k/(2k+1)!.
 # |z|<=4 => after k=24 tail <=2*4^25/50! < 2^-48.
 if abs(z)>4:return None
 zz=enc(z);one=(I,I);C=one;S=one;tc=one;ts=one
 for k in range(1,25):
  tc=idiv(imul(tc,zz),(2*k-1)*(2*k));ts=idiv(imul(ts,zz),(2*k)*(2*k+1));C=iadd(C,tc);S=iadd(S,ts)
 return (C[0]-1,C[1]+1),(S[0]-1,S[1]+1)
def fibre(vs,areas,target=T):
 am=min(areas)
 if am<=0:return False
 H=hull(x for v in vs for x in v)
 if len(H)<3:return False
 # hull integer coords; rescale x and y by global projection length. m transformed by length^2.
 candidates=dirs(vs); rank=[]
 for d in candidates:
  n=norm2(d)
  if not n:continue
  xx=[dot(p,d) for p in H];yy=[cross(d,p) for p in H];w=max(xx)-min(xx);h=max(yy)-min(yy)
  if not w or not h:continue
  val=PI2*n*am/2*(Q(1,w*w)+Q(1,h*h))
  if val>target:return True
  rank.append((val,d,xx,yy,w,n))
 for _,d,xx,yy,w,n in sorted(rank,reverse=True)[:3]:
  verts=[(Q(x-min(xx),w),Q(y-min(yy),w)) for x,y in zip(xx,yy)];U=[];L=[]
  for a,b in zip(verts,verts[1:]+verts[:1]):
   if a[0]==b[0]:continue
   slope=(b[1]-a[1])/(b[0]-a[0]);(U if b[0]<a[0] else L).append((slope,a[1]-slope*a[0]))
  if not U or not L:continue
  m=target*2*w*w/(am*n)*(1+Q(1,10**9))
  for N in (24,48,96):
   h=Q(1,N);phi=(0,0);der=(I,I);ok=True
   for i in range(N):
    a=i*h;b=(i+1)*h
    W=min(max((u[0]-l[0])*a+u[1]-l[1],(u[0]-l[0])*b+u[1]-l[1]) for u in U for l in L)
    v=min(PI2/W**2,m+4*N*N) if W>0 else m+4*N*N;z=(v-m)*h*h
    cs=transfer(z)
    if cs is None:ok=False;break
    C0,S0=cs; phin=iadd(imul(C0,phi),imul(imul(enc(h),S0),der));dern=iadd(imul(imul(enc((v-m)*h),S0),phi),imul(C0,der));phi,der=phin,dern
    if phi[0]<=0:ok=False;break
    # |z|<=4 ensures trig segment length <=2<pi. Nodes positive implies PS.
   if ok:return True
 return False

def check(chart,j,lo,hi,label,table=0):
 D,s,ll,hh,qs=quad(chart,lo,hi);rr=ear_ranges(qs);y,Y=rr[0];h,H=rr[1];k,K=rr[2];d,Dmax=rr[3];t0,u0=ll[4:];t1,u1=hh[4:];tu=Q(t0*u0,D*D)
 if any(rr[i][1]<=0 for i in range(4)):return 'PROVED' if label=='V' else 'PROVED' # empty member set
 if label=='V':return None
 if label=='E2':return 'PROVED' if all(q3[1]*4>=s and 8*inv[2]<(s-q3[0])*s for q2,q3,inv,a in qs) else None
 if label=='E3':return 'PROVED' if all(q2[0]*4>=s and 8*inv[3]<q2[0]*s for q2,q3,inv,a in qs) else None
 if label in ('X','TX','T'):
  gs=untouched(j,rr)
  if label=='X':
   if any(tu>g for g in gs) or stub(j,D,ll,hh,qs) or (chart=='E3' and homface(lo,hi,j,tu)):return 'PROVED'
   return None
  lb=lower_quad(chart,s,qs,rr);amin=min(q[3] for q in qs)
  tau=min(Q(1),Q(max(0,Dmax),max(0,y)+max(0,Dmax))) if j==3 and Dmax>0 else min(Q(1),Q(H,amin)) if j==0 and amin>0 else Q(1)
  if tau<=0:tau=Q(1)
  def yes(L):return (L>T and (Q(t1*u1,D*D)<(1-T/L)/tau if label=='T' else any((1-T/L)/tau>g for g in gs)))
  if yes(lb):return 'PROVED'
  if chart=='std' and table and yes(Q(table)):return 'CONDITIONAL_TABLE'
  # A threshold-targeted independent quadrilateral fibre proof.
  cut=Q(t1*u1,D*D) if label=='T' else min(gs,default=Q(2))
  den=1-cut*tau
  if den>0 and amin>0:
   qvs=[((0,0),(s,0),a,b) for a,b,_,_ in qs];areas=[r[3] for r in qs]
   if fibre(qvs,areas,T/den):return 'PROVED'
  return None
 S,vs=pent(D,s,ll,hh,qs,j)
 if label=='N':return 'PROVED' if geom_N(vs) else None
 areas=[area2(v) for v in vs]
 if label=='S':return 'PROVED' if strip(vs,areas) else None
 if label=='R':return 'PROVED' if tri_excl(vs,areas) else None
 if label=='O':return 'PROVED' if fibre(vs,areas) else None
 return None
