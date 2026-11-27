"""Offline benchmark on frozen data, reservoir samples and EXACT same old one-per-block cases.
Writes actual completed observations; no extrapolation presented as proof.
"""
import argparse,collections,itertools,json,multiprocessing as mp,os,pathlib,random,time
from driver import STAGES,LG,enumerate_stage,sampled,work,limits,canonical,identity,MAX_ROW
P=pathlib.Path
FOCUS={'E2':{'T','TX','O'},'E3':{'T','TX','O'},'std':{'O'}}
def old_samples(data,old):
 # Old runner used per-block, per-label random.Random(20261122).shuffle(v), select first.
 # Only blocks/categories actually reported in B_spot_fixed.jsonl are used.
 roots={c:json.loads((data/'out'/c/'roots.json').read_text())['blocks'] for c in FOCUS}
 for line in (old/'B_spot_fixed.jsonl').open():
  r=json.loads(line);c,b=r['chart'],r['block'];bad={(x['line'],x['label']) for x in r['unproved']};by=collections.defaultdict(list)
  lbs={tuple(k):L for k,lo,hi,L in roots[c][b]}
  for i,l in enumerate((data/'out'/c/'blocks'/b/'leaves.jsonl').open(),1):
   k,j,lo,hi,lab=json.loads(l)
   if lab in FOCUS[c] and lab in r['counts']:by[lab].append((i,(k,j,lo,hi,lab)))
  for lab,rr in by.items():
   assert r['counts'][lab]['checked']==1
   random.Random(20261122).shuffle(rr);i,(k,j,lo,hi,lab)=rr[0]
   from fractions import Fraction as Q
   t=dict(id=f'B/{c}/{b}/{i}',chart=c,block=b,key=k,j=j,lo=lo,hi=hi,label=lab,table=str(Q(lbs[tuple(k)])) if c=='std' else '0',stage='B_'+c)
   old_status='UNPROVED' if (i,lab) in bad else 'PROVED'
   assert (old_status=='UNPROVED')==bool(r['counts'][lab].get('unproved'))
   yield f'focus_{c}_{lab}',t,old_status

def execute(x):
 group,t,old,engine,budget=x;r=work((t,engine,budget));return group,t,old,r

def main():
 p=argparse.ArgumentParser();p.add_argument('data',type=P);p.add_argument('old_logs',type=P);p.add_argument('out',type=P);p.add_argument('--jobs',type=int,default=4);p.add_argument('--seconds',type=int,default=0);p.add_argument('--sample',type=int,default=200);p.add_argument('--qsample',type=int,default=1000);p.add_argument('--budget',type=int,default=64);p.add_argument('--focus-only',action='store_true');p.add_argument('--bench-only',action='store_true');p.add_argument('--groups',nargs='+',help='glob patterns for group names');p.add_argument('--old-fail-first',action='store_true');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 groups=collections.defaultdict(list)
 if not a.focus_only:
  for s in STAGES:
   ts=sampled(enumerate_stage(a.data,s),a.sample)
   groups['bench_'+s].extend((t,None) for t in ts)
 if not a.bench_only:
  for g,t,old in old_samples(a.data,a.old_logs):groups[g].append((t,old))
  analytic={json.loads(l)['id'] for l in (a.old_logs/'q_analytic.jsonl').open()}
  q=sampled((t for t in enumerate_stage(a.data,'q_other') if t['method'] in ('fibre','root_fibre') and t['id'] not in analytic),a.qsample)
  groups['q_fibre1000']=[(t,None) for t in q]
 if a.groups:
  import fnmatch
  groups={g:v for g,v in groups.items() if any(fnmatch.fnmatch(g,p) for p in a.groups)}
 counts={g:collections.Counter() for g in groups};done={};files={}
 for g,v in groups.items():
  # Deterministic random within each group; old failure cases are not preferentially selected.
  random.Random(20261122).shuffle(v)
  if a.old_fail_first:v.sort(key=lambda x:x[1]!='UNPROVED')
  path=a.out/(g+'.jsonl');d={}
  if path.exists():
   for l in path.open():
    r=json.loads(l);assert r['task_sha256']==identity(r['task']);d[r['id']]=r;counts[g][r['status']]+=1
  done[g]=d;files[g]=path.open('a',buffering=1)
 print('PLAN',json.dumps({g:len(v) for g,v in groups.items()},sort_keys=True),flush=True)
 def one(g,v):
  for t,o in v:
   if t['id'] not in done[g]:yield g,t,o,'scipy',a.budget
 # Round robin so every requested class receives work even if runtime capped.
 it=(v for row in itertools.zip_longest(*(one(g,v) for g,v in groups.items())) for v in row if v is not None)
 st=time.monotonic();completed=0
 with mp.get_context('spawn').Pool(a.jobs,initializer=limits,initargs=(1024*1024**2,),maxtasksperchild=128) as pool:
  rr=pool.imap_unordered(execute,it,chunksize=1)
  while True:
   try:g,t,o,r=rr.next(timeout=max(.01,a.seconds-(time.monotonic()-st)) if a.seconds else None)
   except StopIteration:break
   except mp.TimeoutError:print('TIME LIMIT',flush=True);break
   r['task']=t;r['old_status']=o;files[g].write(canonical(r)+'\n');counts[g][r['status']]+=1;completed+=1
   if completed%25==0:print('COMPLETED',completed,'WALL',round(time.monotonic()-st,2),'COUNTS',json.dumps(counts,sort_keys=True),flush=True)
   if a.seconds and time.monotonic()-st>=a.seconds:break
 for f in files.values():f.close()
 print('SUMMARY',json.dumps(counts,sort_keys=True),'WALL',time.monotonic()-st,flush=True)
if __name__=='__main__':main()
