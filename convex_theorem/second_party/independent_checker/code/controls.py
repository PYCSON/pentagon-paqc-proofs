#!/usr/bin/env python3
import json,copy,time,random
from fractions import Fraction as Q
from lg_search import *
from pathlib import Path

def point(name,v,n):
 t=time.monotonic();u,c=proposal('point',[],[],0,n,v);d,p,q,l,h=assemble('point',[],[],0,n,u,c,TARGET,v)
 assert l is not None and l<=h
 print('POINT',name,n,float(l),float(h),'wall',time.monotonic()-t,flush=True)
 return l,h
for name,v in [('square',[(0,0),(1,0),(1,1),(0,1)]),('rational near-equilateral triangle',[(0,0),(1,0),(Q(1,2),Q('0.8660254037844386'))]),('12-gon on approximate circle',[(Q(round(10**8*np.cos(2*np.pi*i/12)),10**8),Q(round(10**8*np.sin(2*np.pi*i/12)),10**8)) for i in range(12)])]:point(name,v,8)
r=random.Random(25)
for i in range(5):
 v=[(Q(round((np.cos(2*np.pi*k/5)+r.uniform(-.07,.07))*10**8),10**8),Q(round((np.sin(2*np.pi*k/5)+r.uniform(-.07,.07))*10**8),10**8)) for k in range(5)];point('random pentagon '+str(i),v,8)
rows=[json.loads(l) for l in open('logs/near_sample.jsonl')];w=next(r['witness'] for r in rows if r['status']=='PROVED');assert replay(w);print('BASELINE',w['chart'],w['N'],w['digest'],flush=True)
bad=copy.deepcopy(w);bad['u'][0][0]+=1;print('COEFFICIENT_TAMPER_REJECT',not replay(bad),flush=True);assert not replay(bad)
lo=list(map(Q,w['lo']));hi=list(map(Q,w['hi']));u=w['u'];chi=w['chi'];ds,pd,pc,lb,ub=assemble(w['chart'],lo,hi,w['j'],w['N'],u,chi,Q(19));assert ub<19;tree,n=treeproof(ds+[pd,pc],6,64);assert tree is None;print('TARGET19_REJECT',float(ub),flush=True)
L=[(3*a-b)/2 for a,b in zip(lo,hi)];H=[(3*b-a)/2 for a,b in zip(lo,hi)];u2=[u[0]]+[[2*x for x in v] for v in u[1:]];c2=[chi[0]]+[[2*x for x in v] for v in chi[1:]];ds,pd,pc,*_=assemble('near',L,H,3,w['N'],u2,c2,TARGET);tree,n=treeproof(ds+[pd,pc],6,64);print('DOUBLED_BOX', 'REJECTED' if tree is None else 'PROVED', 'parts',n,flush=True)
# Polynomial identity cross-route: evaluate affine family at a rational corner, then integrate at that point.
def identity(w):
 lo=list(map(Q,w['lo']));hi=list(map(Q,w['hi']));d=len(lo);z=[1 if i%2==0 else -1 for i in range(d)];pt=[b if s==1 else a for a,b,s in zip(lo,hi,z)];u=w['u'];chi=w['chi'];uu=[[a+sum(s*row[i] for s,row in zip(z,u[1:])) for i,a in enumerate(u[0])]];cc=[[a+sum(s*row[i] for s,row in zip(z,chi[1:])) for i,a in enumerate(chi[0])]]
 ds,pd,pc,*_=assemble(w['chart'],lo,hi,w['j'],w['N'],u,chi,Q(w['L']));gg,S,vs=geometry(w['chart'],lo,hi,w['j']);O=[sum(v[i].get(0,0) for v in vs)/len(vs) for i in range(2)];p1=prod(ev(x,z) for x in ds);ds2,pd2,pc2,*_=assemble(w['chart'],pt,pt,w['j'],w['N'],uu,cc,Q(w['L']),origin=O);_,S2,_=geometry(w['chart'],pt,pt,w['j'],origin=O);p2=prod(x.get(0,0) for x in ds2);assert Q(ev(pd,z),S*S*p1)==Q(pd2.get(0,0),S2*S2*p2);assert Q(ev(pc,z),S**4*p1)==Q(pc2.get(0,0),S2**4*p2);print('IDENTITY',w['chart'],len(pd),len(pc),'PASS',flush=True)
for file in ['near_sample','far_sample','q_sample']:
 r=next(json.loads(l)['witness'] for l in open('logs/'+file+'.jsonl') if json.loads(l)['status']=='PROVED');identity(r)
