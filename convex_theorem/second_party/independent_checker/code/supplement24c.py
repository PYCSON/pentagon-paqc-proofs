#!/usr/bin/env python3
"""Retry ONLY r3 UNPROVED/ERROR; original task identity remains unchanged."""
import argparse, fcntl, json, multiprocessing as mp, os, pathlib, tempfile, time, sys
import driver as D
import replay24c as R
Path=pathlib.Path

def retry(args):
 t,engine,budget=args;u=dict(t)
 if t['stage'] in D.LG:u['N']=min(12,t['N']+2)
 result=D.work((u,engine,budget));result['task_sha256']=D.identity(t)
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,default=D.BASE/'work');p.add_argument('--results',type=Path);p.add_argument('--data-zip',type=Path,default=D.BASE/'PENT24_bundle_b8ab7c1.zip');p.add_argument('--stage',choices=['all','A','B',*D.STAGES],default='all');p.add_argument('--number',type=int,default=1);p.add_argument('--engine',choices=['scipy','stdlib'],default='scipy');p.add_argument('--budget',type=int,default=128);p.add_argument('--jobs',type=int,default=20);p.add_argument('--disk-min',type=float,default=20);p.add_argument('--mem-min',type=float,default=2);p.add_argument('--task-timeout',type=float,default=1800);a=p.parse_args()
 if not 1<=a.jobs<=20 or not 64<=a.budget<=256 or a.number<1 or min(a.disk_min,a.mem_min)<0 or a.task_timeout<=0:raise ValueError('invalid supplement options')
 R.check_code();work=a.work.resolve();root=(a.results or work/'results').resolve();dest=work/f'supplement_{a.number}';work.mkdir(exist_ok=True,parents=True)
 with (work/'run.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);D.guard(work,a);plan=R.load_plan(root)
  if root.parent!=work or root.name!=('sample_'+str(plan['config']['sample']) if plan['config']['sample'] else 'results'):raise ValueError('--work must be the original generation work parent')
  if a.engine=='scipy':
   import importlib.util
   if importlib.util.find_spec('numpy') is None or importlib.util.find_spec('scipy') is None:raise ValueError('SciPy absent; use --engine stdlib')
  cfg=dict(version='PENT24c-supplement-v1',base_code=R.R3,code=R.supplement_id(),main_plan_sha256=D.sha(root/'plan.json'),zip_sha256=D.ZIP_SHA,engine=a.engine,budget=a.budget,degree_add=2)
  dest.mkdir(exist_ok=True);pp=dest/'plan.json';new={'config':cfg,'main_totals':plan['totals'],'sample':plan['config']['sample']}
  if pp.exists() and json.loads(pp.read_text())!=new:raise ValueError('supplement config changed; choose another --number')
  D.atomic(pp,new)
  # Once supplements exist, summary must rerun affected replay stages to reflect them.
  with tempfile.TemporaryDirectory(prefix='supp-data-',dir=work) as td:
   data=D.prepare(a.data_zip,Path(td));D.guard(work,a)
   for s in R.selected(a.stage):
    counts={};n=0;mem=1024 if s in ('far_E2','far_E3','B_E2','B_E3') else 512;jobs=min(a.jobs,10240//mem,10 if mem==1024 else 20)
    with mp.get_context('spawn').Pool(jobs,initializer=D.limits,initargs=(mem*1024**2,),maxtasksperchild=128) as pool:
     for idx,batch in enumerate(D.chunks(D.task_iter(data,s,plan['config']['sample']))):
      D.guard(work,a);expected={t['id']:t for t in batch};rel=Path(s)/f'{idx:06d}.jsonl';mainpath=root/rel
      main=D.load_shard(mainpath,expected,False,True);eligible={i:expected[i] for i,r in main.items() if r['status'] in ('UNPROVED','ERROR')};n+=len(eligible)
      path=dest/rel;done=D.load_shard(path,eligible,True,True)
      for v in done.values():counts[v['status']]=counts.get(v['status'],0)+1
      todo=( (t,a.engine,a.budget) for i,t in eligible.items() if i not in done )
      if len(eligible)==len(done):continue
      path.parent.mkdir(parents=True,exist_ok=True)
      with path.open('ab',buffering=0) as f:
       for r in R.bounded_map(pool,retry,todo,work,a):
        if r['status']=='ERROR' and 'MemoryError' in r.get('reason',''):raise InterruptedError('supplement worker memory limit')
        f.write((D.canonical(r)+'\n').encode());os.fsync(f.fileno());counts[r['status']]=counts.get(r['status'],0)+1
      D.log('SUPPLEMENT SHARD',s,idx,counts)
    D.atomic(dest/(s+'_summary.json'),dict(stage=s,selected=n,counts=counts,updated=time.time()));D.log('SUPPLEMENT STAGE',s,n,counts)
if __name__=='__main__':
 try:main()
 except (InterruptedError,KeyboardInterrupt,BlockingIOError) as ex:D.log('SAFE STOP',repr(ex));sys.exit(75)
