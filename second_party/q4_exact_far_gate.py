#!/usr/bin/env python3
"""Original rational geometry gate for convex five-gons (PENT-04 Q4).
No external code used. Handles boundedness of 3/4 halfplane intersections;
computes supporting polygon exact area, and safe rational sufficient gates.
This is a point proof/spec demonstrator, NOT a full interval box cover.
"""
from fractions import Fraction as F
from itertools import combinations
import json,sys
X=F
def p(a):return (X(str(a[0])),X(str(a[1])))
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def area(V):return sum(cross(V[i],V[(i+1)%len(V)]) for i in range(len(V)))/2
def orient(a,b,x):return cross(sub(b,a),sub(x,a))
def intersect(a,b,c,d):
 e=sub(b,a); f=sub(d,c); den=cross(e,f)
 if not den:return None
 t=cross(sub(c,a),f)/den
 return (a[0]+t*e[0],a[1]+t*e[1])
def convex_hull(points):
 pts=sorted(set(points))
 if len(pts)<3:return []
 def turn(a,b,c):return cross(sub(b,a),sub(c,b))
 lo=[];hi=[]
 for pt in pts:
  while len(lo)>1 and turn(lo[-2],lo[-1],pt)<=0:lo.pop()
  lo.append(pt)
 for pt in pts[::-1]:
  while len(hi)>1 and turn(hi[-2],hi[-1],pt)<=0:hi.pop()
  hi.append(pt)
 return lo[:-1]+hi[:-1]
def recession_exists(edges):
 dirs=[]
 for a,b in edges:
  e=sub(b,a);dirs.extend([e,(-e[0],-e[1])])
 return any(all(cross(sub(b,a),d)>=0 for a,b in edges) for d in dirs)
def outer_polygon(V,sel):
 edges=[(V[i],V[(i+1)%5]) for i in sel]
 if recession_exists(edges):return None
 cand=[]
 for (a,b),(c,d) in combinations(edges,2):
  z=intersect(a,b,c,d)
  if z is not None and all(orient(x,y,z)>=0 for x,y in edges):cand.append(z)
 Q=convex_hull(cand)
 if len(Q)<3 or area(Q)<=0:return None
 assert all(all(orient(a,b,z)>=0 for a,b in edges) for z in Q)
 assert all(all(orient(a,b,z)>=0 for a,b in edges) for z in V)
 return Q

def examine(V):
 A=area(V);L=[sum(z*z for z in sub(V[(i+1)%5],V[i])) for i in range(5)]
 if A<=0:return {'status':'non_ccw'}
 if not all(orient(V[i],V[(i+1)%5],V[j])>0 for i in range(5) for j in range(5) if j not in (i,(i+1)%5)):
  return {'status':'non_strict_convex'}
 if max(L)>1:return {'status':'not_max_edge_normalization'}
 assert L[0]==1
 # Use conservative threshold rational with huge safety margin against producer Fstar_up.
 b4=(F(9586,10000));b3=F(8302,10000);HP=F(3068,100)
 per_ub=0
 # Upper bounds of perimeter are unnecessary, use exact rational lower bound sqrt(L) via bisection;
 # or use (sum exact L_i)^2/5 <= Per^2 (Cauchy lower bound) as cheap sufficient gate.
 # Perron lower rational via integer isqrt scaled to 2^-60.
 from math import isqrt
 floor=[]
 for l in L:
  q=(l.numerator<<120)//l.denominator
  floor.append(F(isqrt(q),1<<60))
 perlo=sum(floor)
 HPpasses=perlo*perlo/A>HP
 candidates=[]
 for n in (3,4):
  for sel in combinations(range(5),n):
   Q=outer_polygon(V,sel)
   if Q is None:continue
   rat=A/area(Q)
   threshold=b3 if n==3 else b4
   candidates.append((n,sel,rat,rat>threshold))
 return {'status':'convex', 'area':str(A), 'per2_over_A_lower':str(perlo*perlo/A),
         'HP_pass':HPpasses,'best3':str(max((z[2] for z in candidates if z[0]==3),default=F(0))),
         'best4':str(max((z[2] for z in candidates if z[0]==4),default=F(0))),
         'halfplane_pass':any(z[3] for z in candidates), 'bounded_subsets':len(candidates)}

def run():
 # Example normalized regular P5 and deformed near-delta rational instances.
 examples=[[(0,0),(1,0),(F(1309,1000),F(951,1000)),(F(1,2),F(1539,1000)),(-F(309,1000),F(951,1000))],
           [(0,0),(1,0),(F(9,10),F(1,10)),(F(7,10),F(1,5)),(0,F(1,5))]]
 for i,raw in enumerate(examples):
  V=[p(z) for z in raw];print('example',i,examine(V))
 print('SOUNDNESS: each accepted 3/4 subset has exactly computed bounded intersection Q; point-only, no box proof.')
if __name__=='__main__':run()

