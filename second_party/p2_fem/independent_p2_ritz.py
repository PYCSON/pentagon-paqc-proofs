#!/usr/bin/env python3
"""Original P2 conforming FEM Rayleigh upper-bound trial on an inscribed pentagon.
No bundled source executed. Numerical eigenvector only chooses trial coefficients;
independent Decimal rational interval evaluator will verify its Rayleigh quotient.
"""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree
import sympy as sy
from pathlib import Path
import json,sys,time

def reference_matrices():
    x,y=sy.symbols('x y'); l=[1-x-y,x,y]
    b=[t*(2*t-1) for t in l]+[4*l[0]*l[1],4*l[1]*l[2],4*l[2]*l[0]]
    def integral(poly):
        p=sy.Poly(sy.expand(poly),x,y);return sum(a*sy.factorial(i)*sy.factorial(j)/sy.factorial(i+j+2) for (i,j),a in p.terms())
    M=np.array([[float(integral(f*g)) for g in b] for f in b]);X=[sy.diff(f,x) for f in b];Y=[sy.diff(f,y) for f in b]
    xx=np.array([[float(integral(f*g)) for g in X] for f in X]);xy=np.array([[float(integral(f*g)) for g in Y] for f in X]);yy=np.array([[float(integral(f*g)) for g in Y] for f in Y]);return M,xx,xy,yy

def mesh(m,s=1-1e-10):
    angle=2*np.pi*np.arange(5)/5; V=np.stack([s*np.cos(angle),s*np.sin(angle)],1)
    nodes=[]; idx={};tr=[]
    def put(p):
        key=tuple(np.round(p,14));
        if key not in idx:idx[key]=len(nodes);nodes.append(np.array(p))
        return idx[key]
    # barycentric wedge grid
    for k in range(5):
        idxs={}
        for a in range(m+1):
            for b in range(m+1-a):
                idxs[(a,b)]=put((a*V[k]+b*V[(k+1)%5])/m)
        for a in range(m):
            for b in range(m-a):
                tr.append([idxs[(a,b)],idxs[(a+1,b)],idxs[(a,b+1)]])
                if a+b<m-1:tr.append([idxs[(a+1,b)],idxs[(a+1,b+1)],idxs[(a,b+1)]])
    nodes=np.array(nodes);tr=np.array(tr,dtype=np.int64)
    edges={};mid=[];els=[]
    for tri in tr:
        mids=[]
        for a,b in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]:
            key=tuple(sorted((int(a),int(b))))
            if key not in edges:edges[key]=len(nodes)+len(mid);mid.append((nodes[a]+nodes[b])/2)
            mids.append(edges[key])
        els.append([*tri,*mids])
    pts=np.vstack([nodes,np.array(mid)]);els=np.array(els,dtype=np.int64)
    # boundary point test: point lies on any of five outer lines
    boundary=np.zeros(len(pts),bool)
    for k in range(5):
        a,b=V[k],V[(k+1)%5]; e=b-a;cr=e[0]*(pts[:,1]-a[1])-e[1]*(pts[:,0]-a[0]);t=((pts-a)@e)/(e@e)
        boundary |= (abs(cr)<1e-10)&(t>=-1e-10)&(t<=1+1e-10)
    return V,pts,els,boundary

def assemble(V,pts,els,boundary):
    M,xx,xy,yy=reference_matrices(); inds=np.where(~boundary)[0]; idmap=np.full(len(pts),-1,dtype=int);idmap[inds]=np.arange(len(inds))
    ri=[];ci=[];ka=[];ma=[]
    for t in els:
        a,b,c=pts[t[:3]];B=np.stack([b-a,c-a],1);det=abs(np.linalg.det(B)); Bi=np.linalg.inv(B)
        g=Bi@Bi.T
        K=det*(g[0,0]*xx+g[0,1]*(xy+xy.T)+g[1,1]*yy)
        MM=det*M
        for i in range(6):
            ii=idmap[t[i]]
            if ii<0:continue
            for j in range(6):
                jj=idmap[t[j]]
                if jj<0:continue
                ri.append(ii);ci.append(jj);ka.append(K[i,j]);ma.append(MM[i,j])
    n=len(inds)
    K=sp.coo_matrix((ka,(ri,ci)),shape=(n,n)).tocsr()
    M=sp.coo_matrix((ma,(ri,ci)),shape=(n,n)).tocsr()
    return K,M,inds

def run(m,path):
    t0=time.time();V,pts,els,boundary=mesh(m);K,M,inds=assemble(V,pts,els,boundary)
    vals,vecs=eigsh(K,M=M,k=1,sigma=7.957,which='LM',tol=1e-11)
    v=vecs[:,0];q=float(v@K@v/(v@M@v));area=0.5*sum(np.linalg.det(np.stack([V[i],V[(i+1)%5]],1)) for i in range(5));F=q*area
    print('m',m,'points',len(pts),'internal',len(inds),'elements',len(els),'eigapprox',q,'Fapprox',F,'time',round(time.time()-t0,1),flush=True)
    np.savez_compressed(path,V=V,pts=pts,els=els,boundary=boundary,trial=v,ids=inds)
    return F
if __name__=='__main__':
    m=int(sys.argv[1]);run(m,sys.argv[2]);

