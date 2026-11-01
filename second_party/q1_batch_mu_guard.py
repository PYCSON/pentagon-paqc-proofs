#!/usr/bin/env python3
"""Independent checker driver, resumable: does not execute bundled producer source.
Input MPS parameters are treated strictly as data. A PASS requires 5 exhaustive edges.
"""
import argparse,concurrent.futures,hashlib,json,pathlib,re,subprocess,sys,time,math
from mpmath import iv
ROOT=pathlib.Path(__file__).resolve().parent
FOLD=ROOT/'lines/D-topics/probes/PENT/G1/replay_corners/all'
FILES=['mps_cell_liu_point_rho0.025_v6.json','mps_cell_r0.1_id71_rho0.012_v6.json','mps_t1_r0.1_id71_rho0.015_rigorous_v6.json','mps_t1_liu_point_rho0.03_rigorous_v6.json','mps_t1_id301_rho0.0239_rigorous_v6.json']
P=ROOT/'q1_checkpoints';P.mkdir(exist_ok=True)
def job(item):
 f,c,e,n=item
 key=f'{f[:-5]}_corner{c:02d}_edge{e}_N{n}'
 target=P/(key+'.json');inp=P/(f'{f[:-5]}_corner{c:02d}_input.json')
 if target.exists():return key,'CACHED'
 cmd=[sys.executable,str(ROOT/'independent_mps_iv.py'),str(inp),str(n),str(e)]
 try:
  run=subprocess.run(cmd,text=True,capture_output=True,timeout=7200)
  (P/(key+'.log')).write_text('COMMAND '+repr(cmd)+'\nSTDOUT\n'+run.stdout+'\nSTDERR\n'+run.stderr)
  m=re.search(r'edge\s+'+str(e)+r'\s+sup\s+([\d.eE+\-]+)',run.stdout)
  if run.returncode or not m: return key,'FAIL:'+str(run.returncode)
  target.write_text(json.dumps({'sup':m.group(1),'log':key+'.log','N':n,'edge':e}))
  return key,'OK'
 except subprocess.TimeoutExpired as ex:
  (P/(key+'.timeout')).write_text(str(ex))
  return key,'TIMEOUT'
def summary(f,c,n):
 data=json.loads((P/(f'{f[:-5]}_corner{c:02d}_input.json')).read_text())
 parts=[P/(f'{f[:-5]}_corner{c:02d}_edge{e}_N{n}.json') for e in range(5)]
 if not all(p.exists() for p in parts):return {'file':f,'corner':c,'status':'UNRESOLVED','edges':sum(p.exists() for p in parts)}
 import independent_mps_iv as M
 sd=M.setup(data);V,k,coef,fr=sd;x0=[sum(x[j] for x in V)/5 for j in (0,1)]
 radii=[]
 for j in range(5):
  a=V[j];b=V[(j+1)%5];e=M.diff(b,a);r=M.cross(e,M.diff(x0,a))/iv.sqrt(M.dot(e,e));radii.append(M.low(r))
 R=iv.mpf(str(min(radii)*0.99));u0=M.evaluate(x0,None,sd)[0];J1=M.bessel_series(iv.mpf(1),k*R,iv.mpf(1))[0]
 l2=M.low((u0*2*M.P*R*J1/k)**2/(M.P*R*R))
 sup=max(float(json.loads(p.read_text())['sup']) for p in parts)
 area=iv.mpf(0)
 for j in range(5): area+=M.cross(V[j],V[(j+1)%5])/2
 mu=iv.sqrt(area)*iv.mpf([str(sup),str(math.nextafter(sup,math.inf))])/iv.sqrt(iv.mpf(str(l2)))
 lam=M.Q(data['lam_tilde']); flo=M.low(lam/(1+mu)*area)
 mu_hi=M.high(mu)
 # FHM gives a finite upper eigenvalue bound only for mu < 1.
 # Fail closed before ever dividing by 1-mu.
 mu_ok=(mu_hi < 1.0)
 lhs=M.high(lam/(1-mu)*area) if mu_ok else None
 producer=float(data['producer_F_lo_v6']);margin=flo-producer
 hks_ok=bool(mu_ok and lhs < 36.33682907107446)
 status='PASS' if margin>=0 and hks_ok else 'FAIL'
 return {'file':f,'corner':c,'status':status,'F_lo':flo,'producer_F_lo_v6':producer,'margin':margin,'mu_upper':mu_hi,'mu_lt_1':mu_ok,'HKS_left_upper':lhs,'HKS_pass':hks_ok,'N':n,'sup':sup,'L2_lower':l2}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=300);ap.add_argument('--corners',type=int,default=64);ap.add_argument('--files',type=int,default=2);ap.add_argument('--workers',type=int,default=3);ap.add_argument('--edge-limit',type=int,default=5);args=ap.parse_args()
 items=[]
 for f in FILES[:args.files]:
  data=json.loads((FOLD/f).read_text());assert len(data['corners'])==64
  for c in range(min(args.corners,64)):
   inp=P/(f'{f[:-5]}_corner{c:02d}_input.json');inp.write_text(json.dumps(data['corners'][c]))
   for e in range(args.edge_limit):items.append((f,c,e,args.N))
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
  for i,(key,result) in enumerate(pool.map(job,items),1):print(i,len(items),key,result,flush=True)
 results=[summary(f,c,args.N) for f in FILES[:args.files] for c in range(min(args.corners,64))]
 (P/f'summary_N{args.N}.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
