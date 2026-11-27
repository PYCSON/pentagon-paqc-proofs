#!/usr/bin/env python3
"""Only 13 named near targets, after full r3 and first normal supplement. No in-place result changes."""
import argparse,fcntl,json,multiprocessing as mp,os,sys,tempfile
from pathlib import Path
import driver as D
import replay26a as R
import outside26a as O

def job(args):return O.search(*args)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,default=D.BASE/'work');p.add_argument('--data-zip',type=Path,default=D.BASE/'PENT24_bundle_b8ab7c1.zip');p.add_argument('--number',type=int,default=2);p.add_argument('--engine',choices=('scipy','stdlib'),default='scipy');p.add_argument('--jobs',type=int,default=4);p.add_argument('--disk-min',type=float,default=20);p.add_argument('--mem-min',type=float,default=2);p.add_argument('--task-timeout',type=float,default=7200);a=p.parse_args()
 if not 1<=a.jobs<=20 or a.number<2 or min(a.disk_min,a.mem_min)<0 or a.task_timeout<=0:raise ValueError('bad options')
 work=a.work.resolve();root=work/'results';dest=work/f'supplement_{a.number}';R.check_code()
 with (work/'run.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);D.guard(work,a);plan=R.load_plan(root)
  if plan['config']['sample']:raise ValueError('not full r3')
  prior=R.supplement_plans(work,root)
  if 'supplement_1' not in prior or prior['supplement_1'][1]['version']!='PENT24c-supplement-v1' or prior['supplement_1'][1]['budget']!=128:raise ValueError('run first normal supplement budget 128 first')
  cfg={'version':'PENT26a-outside-v1','base_code':R.R3,'code':R.outside_id(),'main_plan_sha256':D.sha(root/'plan.json'),'zip_sha256':D.ZIP_SHA,'engine':a.engine,'budget':256,'degree_add':2,'ids':sorted(O.IDS)};new={'config':cfg,'main_totals':plan['totals'],'sample':0}
  # Freeze source eligibility by scanning each full main shard against authenticated task enumeration. No reliance on summary counts.
  with tempfile.TemporaryDirectory(prefix='outside-data-',dir=work) as td:
   data=D.prepare(a.data_zip,Path(td));targets={};found=set()
   for s in D.STAGES:
    total=0
    for idx,batch in enumerate(D.chunks(D.task_iter(data,s,0))):
     D.guard(work,a);exp={t['id']:t for t in batch};rel=Path(s)/f'{idx:06d}.jsonl';main=D.load_shard(root/rel,exp,False,True)
     if set(main)!=set(exp):raise ValueError('main incomplete '+str(rel))
     total+=len(main)
     eligible={i:exp[i] for i,r in main.items() if r['status'] in ('UNPROVED','ERROR')}
     first=D.load_shard(prior['supplement_1'][0]/rel,eligible,False,True)
     if set(first)!=set(eligible):raise ValueError('first normal supplement incomplete '+str(rel))
     if s!='near':continue
     sel=O.IDS.intersection(exp);found.update(sel)
     for i in sel.intersection(eligible):
      if first[i]['status'] in ('UNPROVED','ERROR'):targets.setdefault(rel,{})[i]=exp[i]
    if total!=R.TOTALS[s]:raise ValueError('main count mismatch '+s)
   if found!=O.IDS:raise ValueError('13 IDs not matched')
   # Prior rounds may already have proved a target. Avoid unnecessary retries; they will be replayed separately.
   for rel,eligible in list(targets.items()):
    for k,(path,c) in prior.items():
     if k in ('supplement_1',dest.name):continue
     if (path/rel).exists():
      for line in (path/rel).open():
       r=json.loads(line);i=r['id']
       if i in eligible and r['status']=='PROVED' and r['task_sha256']==D.identity(eligible[i]):del eligible[i]
   dest.mkdir(exist_ok=True);pp=dest/'plan.json'
   if pp.exists() and json.loads(pp.read_text())!=new:raise ValueError('supplement plan changed')
   D.atomic(pp,new);D.log('OUTSIDE SELECTED',sum(map(len,targets.values())))
   # Use 1 GiB per worker for N=10; at most 10 GiB. Existing global guard retained.
   with mp.get_context('spawn').Pool(min(a.jobs,10),initializer=D.limits,initargs=(1024**3,),maxtasksperchild=1) as pool:
    for rel,eligible in targets.items():
     path=dest/rel;done=D.load_shard(path,eligible,True,True);todo=((t,a.engine,256) for i,t in eligible.items() if i not in done)
     if len(done)==len(eligible):continue
     path.parent.mkdir(exist_ok=True,parents=True)
     with path.open('ab',buffering=0) as f:
      for r in R.bounded_map(pool,job,todo,work,a):
       if r['status']=='ERROR' and 'MemoryError' in r.get('reason',''):raise InterruptedError('worker memory cap')
       f.write((D.canonical(r)+'\n').encode());os.fsync(f.fileno());D.log('OUTSIDE RESULT',r['id'],r['status'],r['parts'],r['outside_leaves'],r['seconds'])
   D.log('OUTSIDE DONE',dest)
if __name__=='__main__':
 try:main()
 except (InterruptedError,KeyboardInterrupt,BlockingIOError) as e:D.log('SAFE STOP',repr(e));sys.exit(75)
