#!/usr/bin/env python3
"""Independent FHM aggregator from segment certificates with checked exact coverage.
All edge supremum records are in strict own-script logs; arithmetic directed via iv.
"""
from pathlib import Path
from mpmath import iv
import json,math,re,sys
from independent_mps_iv import setup,evaluate,Q,low,high,cross,diff,dot,bessel_series
iv.dps=45
EVID=Path(__file__).resolve().parents[1]; ROOT=Path('/mnt/data/pent04_workspace/lines/D-topics/probes/PENT/G1/replay_corners')

def run(name,N):
    data=json.loads((ROOT/name).read_text())
    rows=[]
    for edge in range(5):
        spans=[]
        for st in range(0,N,100):
            ed=min(st+100,N)
            logfile=EVID/'logs'/f'{name[:-5]}_edge{edge}_{st:04d}-{ed:04d}_N{N}.log'
            if not logfile.exists():raise ValueError(f'MISSING {logfile.name}')
            s=logfile.read_text()
            m=re.search(r'^edge\s+(\d+)\s+sup\s+(\S+)',s,re.MULTILINE)
            if not m or int(m.group(1))!=edge:raise ValueError(f'INVALID log {logfile.name}')
            spans.append((st,ed,float(m.group(2))))
        assert spans[0][0]==0 and spans[-1][1]==N and all(spans[i][1]==spans[i+1][0] for i in range(len(spans)-1))
        rows.append(math.nextafter(max(x[2] for x in spans)*(1+1e-12),math.inf))
    sup=math.nextafter(max(rows),math.inf)
    V,k,coeff,fr=setup(data);cen=[sum(x[j] for x in V)/5 for j in (0,1)]
    radii=[]
    for j in range(5):
        e=diff(V[(j+1)%5],V[j]);r=cross(e,diff(cen,V[j]))/iv.sqrt(dot(e,e))
        assert low(r)>0
        radii.append(low(r))
    R=iv.mpf(str(min(radii)*0.99))
    u0=evaluate(cen,None,(V,k,coeff,fr))[0]
    J=bessel_series(iv.mpf(1),k*R,iv.mpf(1),N=38)[0]
    l2=(u0*2*iv.pi*R*J/k)**2/(iv.pi*R*R)
    l2_lower=low(l2)
    area=sum(cross(V[j],V[(j+1)%5]) for j in range(5))/2
    mu=iv.sqrt(area)*iv.mpf(str(sup))/iv.sqrt(iv.mpf(str(l2_lower)))
    val=Q(data['lam_tilde'])/(1+mu)*area
    fdown=low(val)
    index_up=high(Q(data['lam_tilde'])/(1-mu)*area)
    out=dict(name=name,N=N,edge_sup=rows,sup_upper=sup,l2_lower=l2_lower,mu_upper=high(mu),F_lower=fdown,producer_F_lo=data['producer_F_lo'],margin=fdown-float(data['producer_F_lo']),index_up=index_up,index_ok=index_up<36.33682907107446,accept=fdown>=float(data['producer_F_lo']))
    print(json.dumps(out,indent=2));(EVID/'logs'/f'mps_aggregate_{name[:-5]}_N{N}.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]))

