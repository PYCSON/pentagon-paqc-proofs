#!/usr/bin/env python3
"""PENT24b streaming, resumable local runner. Linux/WSL2; no network, no producer code."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import argparse,collections,fcntl,hashlib,itertools,json,math,multiprocessing as mp,pathlib,random,resource,shutil,signal,sys,tarfile,tempfile,time,zipfile
from fractions import Fraction as Q
from lg_runner import near_tasks,far_tasks,q_tasks
from q_fibre_runner import tasks as q_other_tasks
from proofs24b import refit_lg,lg_bound,discrete_proof,discrete_bound
P=pathlib.Path
BASE=P(__file__).resolve().parent.parent
ZIP_SHA='1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958'
STAGES=('near','q_lg','q_other','far_std','far_E2','far_E3','B_std','B_E2','B_E3')
LG=set(STAGES[:2]+STAGES[3:6]);CHUNK=128;MAX_ROW=4*1024*1024;GIB=1024**3
# Frozen counts including ALL closed-form q parts, not just the fibre subset.
TOTALS={'near':38280,'q_lg':66138,'q_other':None,'far_std':None,'far_E2':None,'far_E3':None,'B_std':None,'B_E2':None,'B_E3':None}
VERSION='PENT24b-20261122-r3'

def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def canonical(o):return json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False)
def identity(t):return hashlib.sha256(canonical(t).encode()).hexdigest()
def atomic(p,obj):
 p=P(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp')
 with tmp.open('w') as f:f.write(canonical(obj)+'\n');f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def log(*a):print(time.strftime('%Y-%m-%dT%H:%M:%S%z'),*a,flush=True)
def code_id():
 h=hashlib.sha256()
 for name in ('driver.py','proofs24b.py','lg_runner.py','lg_exact.py','lg_search.py','lg_stdlib.py','label_exact.py','q_fibre_runner.py'):
  p=BASE/'code'/name;h.update(p.name.encode()+b'\0'+bytes.fromhex(sha(p)))
 return h.hexdigest()

def data_hashes(data):
 return {str(p.relative_to(data)):sha(p) for p in sorted(data.rglob('*')) if p.is_file() and p.name!='DATA-MANIFEST.json'}
def prepare(zpath,work):
 zpath=P(zpath).resolve();data=P(work)/'data';meta=data/'DATA-MANIFEST.json'
 if sha(zpath)!=ZIP_SHA:raise ValueError('FROZEN ZIP SHA256 mismatch')
 if meta.exists():
  exp=json.loads(meta.read_text());got=data_hashes(data)
  if exp!=got:raise ValueError('extracted data hash mismatch; remove data directory and reprepare')
  return data
 if data.exists():raise ValueError('Incomplete data directory; remove it before retry')
 tmp=P(tempfile.mkdtemp(prefix='extract-',dir=work))
 try:
  with zipfile.ZipFile(zpath) as z:
   for n in ('near/blocks.tgz','far/pent_far5_cert.tar.xz'):
    name='PENT24_bundle_b8ab7c1/data/'+n;archive=tmp/P(n).name
    with z.open(name) as f,archive.open('wb') as o:shutil.copyfileobj(f,o)
    with tarfile.open(archive,'r:*') as tf:
     for member in tf:
      path=P(member.name)
      if path.is_absolute() or '..' in path.parts or not(member.isfile() or member.isdir()):raise ValueError('unsafe tar member')
      tf.extract(member,tmp,filter='data')
    archive.unlink()
  if not (tmp/'blocks').is_dir() or not (tmp/'out/qtable/cells.jsonl').exists():raise ValueError('frozen data layout')
  # No production executable is needed or run.
  atomic(tmp/'DATA-MANIFEST.json',data_hashes(tmp));os.replace(tmp,data)
 except BaseException:shutil.rmtree(tmp,ignore_errors=True);raise
 return data

def label_tasks(far,chart):
 roots=json.loads((far/chart/'roots.json').read_text())['blocks']
 for block,rs in sorted(roots.items()):
  lbs={tuple(k):L for k,lo,hi,L in rs}
  with (far/chart/'blocks'/block/'leaves.jsonl').open() as f:
   for i,line in enumerate(f,1):
    k,j,lo,hi,lab=json.loads(line)
    yield dict(id=f'B/{chart}/{block}/{i}',chart=chart,block=block,key=k,j=j,lo=lo,hi=hi,label=lab,table=str(Q(lbs[tuple(k)])) if chart=='std' else '0')

def enumerate_stage(data,stage):
 far=data/'out'
 if stage=='near':it=near_tasks(data)
 elif stage=='q_lg':it=(t for t in q_tasks(far) if t['method']=='lg')
 elif stage=='q_other':it=q_other_tasks(far)
 elif stage.startswith('far_'):
  chart=stage[4:]
  def fs():
   for p in sorted((far/chart/'blocks').glob('*/cert.jsonl')):
    for i,line in enumerate(p.open(),1):
     q=json.loads(line);from lg_exact import TARGET
     yield dict(id=f'{chart}/{p.parent.name}/{i}',chart=chart,lo=q['lo'],hi=q['hi'],j=q['j'],N=q['N'],L=str(TARGET),block=p.parent.name,key=q['key'])
  it=fs()
 elif stage.startswith('B_'):it=label_tasks(far,stage[2:])
 else:raise ValueError(stage)
 for t in it:
  t['stage']=stage;yield t

def sampled(it,k,seed=20261122):
 # Reservoir sampling gives a deterministic uniform sample without retaining full dataset.
 r=random.Random(seed);v=[]
 for i,t in enumerate(it):
  if i<k:v.append(t)
  else:
   j=r.randrange(i+1)
   if j<k:v[j]=t
 return sorted(v,key=lambda t:t['id'])

def task_iter(data,stage,sample):
 it=enumerate_stage(data,stage);return iter(sampled(it,sample)) if sample else it

def chunks(it):
 while batch:=list(itertools.islice(it,CHUNK)):yield batch

def limits(mem):
 # WSL Linux per-process address-space cap; no child may exceed 512MiB by default.
 resource.setrlimit(resource.RLIMIT_AS,(mem,mem));signal.signal(signal.SIGINT,signal.SIG_IGN)
def work(args):
 t,engine,budget=args;st=time.monotonic();r={'id':t['id'],'task_sha256':identity(t),'status':'UNPROVED'}
 try:
  if t['stage'] in LG:
   w,info=refit_lg(t,engine,budget);r.update(info);r['status']='PROVED' if w else 'UNPROVED'
  else:
   w,info=discrete_proof(t,Q(t.get('table','0')),budget);r.update(info)
  if w:r['witness']=w
 except MemoryError:r={'id':t['id'],'task_sha256':identity(t),'status':'ERROR','reason':'MemoryError: worker memory limit'}
 except Exception as e:r.update(status='ERROR',reason=repr(e)[:400])
 r['seconds']=time.monotonic()-st;r['rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 # Cap per-record, hence a 128-record shard cannot exceed 512MiB.
 if len(canonical(r).encode())+1>MAX_ROW:r={k:r[k] for k in ('id','task_sha256','seconds','rss_kib')};r.update(status='UNPROVED',reason='witness exceeds row cap')
 return r

def memavail():
 m=int(next(l.split()[1] for l in open('/proc/meminfo') if l.startswith('MemAvailable:')))*1024
 for root in (P('/sys/fs/cgroup'),):
  try:
   lim=(root/'memory.max').read_text().strip();cur=int((root/'memory.current').read_text());m=min(m,max(0,int(lim)-cur)) if lim!='max' else m
  except (OSError,ValueError):pass
 return m

def guard(work,a):
 if (work/'STOP').exists():raise InterruptedError('STOP requested')
 if shutil.disk_usage(work).free<a.disk_min*GIB:raise InterruptedError('disk free below configured threshold (default 20 GiB)')
 if memavail()<a.mem_min*GIB:raise InterruptedError('available memory below configured threshold (default 2 GiB)')

def load_shard(path,expected,repair=False,metadata_only=False):
 rows={}
 if not path.exists():return rows
 if path.stat().st_size>CHUNK*MAX_ROW:raise ValueError('oversize shard '+str(path))
 with path.open('rb+' if repair else 'rb') as f:
  while True:
   off=f.tell();line=f.readline(MAX_ROW+1)
   if not line:break
   if len(line)>MAX_ROW:raise ValueError('oversize row')
   if not line.endswith(b'\n'):
    if repair:f.truncate(off);f.flush();os.fsync(f.fileno());log('RECOVER truncated tail',path,off);break
    raise ValueError('truncated row')
   r=json.loads(line);i=r['id']
   if i not in expected or r['task_sha256']!=identity(expected[i]) or i in rows:raise ValueError('duplicate or input binding mismatch '+i)
   if r['status'] not in ('PROVED','CONDITIONAL_TABLE','UNPROVED','ERROR'):raise ValueError('invalid status')
   rows[i]={"status":r["status"]} if metadata_only else r
 return rows

def state_write(work,stage,ct,total,pool=None,note=''):
 pids=[p.pid for p in pool._pool if p.is_alive()] if pool else []
 atomic(work/'progress'/f'{stage}.json',dict(stage=stage,total=total,counts=dict(ct),workers=pids,parent=os.getpid(),updated=time.time(),note=note))

def run_stage(data,work,results,stage,total,a):
 root=results/stage;root.mkdir(parents=True,exist_ok=True);ct=collections.Counter();begin=time.monotonic();log('STAGE',stage,'total',total)
 pool=None;ended=False;seen=0;last=0
 try:
  mem=a.edge_worker_mib if stage in ('far_E2','far_E3') else a.worker_mib
  jobs=min(a.jobs,10240//mem)
  log('RESOURCE',stage,'jobs',jobs,'worker_address_space_MiB',mem)
  guard(work,a);pool=mp.get_context("spawn").Pool(jobs,initializer=limits,initargs=(mem*1024**2,),maxtasksperchild=128)
  for idx,batch in enumerate(chunks(task_iter(data,stage,a.sample))):
   seen+=len(batch);expected={t['id']:t for t in batch};path=root/f'{idx:06d}.jsonl';done=load_shard(path,expected,True,True)
   ct.update(r['status'] for r in done.values());todo=iter(t for t in batch if t['id'] not in done);pending={}
   def submit():
    while len(pending)<jobs*2:
     try:t=next(todo)
     except StopIteration:return
     pending[t['id']]=(pool.apply_async(work_func,((t,a.engine,a.budget),)),time.monotonic())
   guard(work,a);submit()
   with path.open('ab',buffering=0) as f:
    while pending:
     guard(work,a);now=time.monotonic();ready=[i for i,(v,s) in pending.items() if v.ready()]
     if not ready:
      if any(now-s>a.task_timeout for v,s in pending.values()):raise InterruptedError('task timeout or worker died; resume with larger timeout or smaller budget in NEW run')
      if now-last>5:state_write(work,stage,ct,total,pool);last=now
      time.sleep(.1);continue
     for i in ready:
      r=pending.pop(i)[0].get()
      if r['status']=='ERROR' and ('MemoryError' in r.get('reason','') or 'failed to map segment' in r.get('reason','')):raise InterruptedError(r['reason'])
      line=(canonical(r)+'\n').encode();f.write(line);ct[r['status']]+=1
      if sum(ct.values())%16==0:os.fsync(f.fileno())
      if time.monotonic()-last>5:state_write(work,stage,ct,total,pool);last=time.monotonic()
     submit()
    os.fsync(f.fileno())
   state_write(work,stage,ct,total,pool)
   if idx%20==0:log(stage,sum(ct.values()),'/',total,dict(ct))
  if seen!=total or sum(ct.values())!=total:raise ValueError(f'coverage count {seen}/{total} rows={sum(ct.values())}')
  ended=True;log('STAGE DONE',stage,dict(ct),'wall',round(time.monotonic()-begin,3))
 finally:
  if pool:pool.terminate();pool.join()
  state_write(work,stage,ct,total,note='DONE' if ended else 'STOPPED; resume')
# Avoid collision with path variable work in run_stage.
work_func=work

def pack(results,work,a):
 guard(work,a);m=results/'SHA256SUMS';m.unlink(missing_ok=True)
 paths=sorted(p for p in results.rglob('*') if p.is_file() and not p.name.endswith('.tmp'))
 size=sum(p.stat().st_size for p in paths)
 # Pessimistic bound: archive no larger than raw size plus tar overhead; allow 20 GiB reserve.
 if shutil.disk_usage(work).free < a.disk_min*GIB+size+len(paths)*1024+10*1024**2:raise InterruptedError('not enough disk space for archive plus 20 GiB reserve')
 with m.open('w') as f:
  for p in paths:f.write(f'{sha(p)}  {p.relative_to(results)}\n')
 dest=work/(('PENT24b_smoke_results.tgz' if a.sample==5 else f'PENT24b_sample_{a.sample}.tgz') if a.sample else 'PENT24b_partial_results.tgz' if a.partial_pack else 'PENT24b_results.tgz');tmp=dest.with_suffix('.tmp')
 try:
  with tarfile.open(tmp,'w:gz',compresslevel=3) as tf:
   for p in paths+[m]:
    guard(work,a);tf.add(p,arcname=str(P('results')/p.relative_to(results)),recursive=False)
  os.replace(tmp,dest);log('ARCHIVE',dest,'SHA256',sha(dest),'bytes',dest.stat().st_size)
 except BaseException:tmp.unlink(missing_ok=True);raise

def run(a):
 if not (1<=a.jobs<=20 and 64<=a.worker_mib<=2048 and 64<=a.edge_worker_mib<=2048 and 1<=a.budget<=64):raise ValueError('jobs 1..20, worker-mib/edge-worker-mib 64..2048, budget 1..64')
 # stage worker allocation <= 10240 MiB + parent <= 10.75 GiB ceiling for this program's address spaces.
 resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,2048*1024**2))
 work=P(a.work).resolve();work.mkdir(parents=True,exist_ok=True)
 with (work/'run.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);(work/'STOP').unlink(missing_ok=True)
  guard(work,a);data=prepare(a.data_zip,work)
  if a.engine=='scipy':
   import importlib.util
   if importlib.util.find_spec('numpy') is None or importlib.util.find_spec('scipy') is None:raise ValueError('SciPy engine absent: install dependencies or use --engine stdlib')
  stages=list(STAGES if a.stage=='all' else [a.stage] if a.stage in STAGES else STAGES[:6] if a.stage=='A' else STAGES[6:]);results=work/('sample_'+str(a.sample) if a.sample else 'results');results.mkdir(exist_ok=True)
  config=dict(version=VERSION,code=code_id(),zip_sha256=ZIP_SHA,sample=a.sample,engine=a.engine,budget=a.budget,chunk=CHUNK)
  plan=results/'plan.json'
  if plan.exists():
   p=json.loads(plan.read_text());assert p['config']==config,'configuration changed: use a fresh work directory'
   totals=p['totals']
  else:
   log('ENUMERATING frozen inventory');full={s:sum(1 for _ in enumerate_stage(data,s)) for s in STAGES};assert all(full[k]==v for k,v in TOTALS.items() if v is not None);assert sum(full[s] for s in STAGES[3:6])==377515;assert sum(full[s] for s in STAGES[6:])==1346929
   totals={s:min(n,a.sample) if a.sample else n for s,n in full.items()};atomic(plan,dict(config=config,totals=totals,full_totals=full));log('INVENTORY',full)
  for s in stages:run_stage(data,work,results,s,totals[s],a)
  summary={}
  for s in STAGES:
   f=work/'progress'/f'{s}.json';summary[s]=json.loads(f.read_text()) if f.exists() else {'stage':s,'total':totals[s],'counts':{},'note':'NOT RUN'}
  atomic(results/'summary.json',summary)
  if all(sum(summary[s]['counts'].values())==totals[s] for s in STAGES):a.partial_pack=False;pack(results,work,a)
  elif a.partial_pack:pack(results,work,a)
  else:log('PARTIAL RUN; other stages not complete; no full results archive')

def progress(a):
 work=P(a.work);p=work/('sample_'+str(a.sample) if a.sample else 'results')/'plan.json';tot=json.loads(p.read_text())['totals'] if p.exists() else TOTALS
 print('stage | finished / total | % | PROVED | CONDITIONAL | UNPROVED/ERROR | active | updated')
 for s in STAGES:
  f=work/'progress'/f'{s}.json';r=json.loads(f.read_text()) if f.exists() else {};c=r.get('counts',{});n=sum(c.values());total=tot[s];pids=r.get('workers',[]);active=0
  for pid in pids:
   try:stat=(P('/proc')/str(pid)/'stat').read_text().split(') ',1)[1].split();active+=int(stat[1])==r.get('parent') and stat[0]!='Z'
   except (OSError,ValueError,IndexError):pass
  print(s,n,'/',total, f'{100*n/total:.2f}%' if total else '?',c.get('PROVED',0),c.get('CONDITIONAL_TABLE',0),c.get('UNPROVED',0)+c.get('ERROR',0),active,time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(r['updated'])) if r else '-',r.get('note',''))
 used=sum(p.stat().st_size for p in work.rglob('*') if p.is_file());print(f'work files {used/GIB:.3f} GiB; disk free {shutil.disk_usage(work if work.exists() else work.parent).free/GIB:.2f} GiB; available memory {memavail()/GIB:.2f} GiB')

def args():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['run','progress']);p.add_argument('--work',type=P,default=BASE/'work');p.add_argument('--data-zip',type=P,default=BASE/'PENT24_bundle_b8ab7c1.zip');p.add_argument('--stage',choices=['all','A','B',*STAGES],default='all');p.add_argument('--jobs',type=int,default=20);p.add_argument('--sample',type=int,default=0);p.add_argument('--engine',choices=['scipy','stdlib'],default='scipy');p.add_argument('--budget',type=int,default=64);p.add_argument('--worker-mib',type=int,default=512);p.add_argument('--edge-worker-mib',type=int,default=1024);p.add_argument('--disk-min',type=float,default=20);p.add_argument('--mem-min',type=float,default=2);p.add_argument('--task-timeout',type=float,default=900);p.add_argument('--partial-pack',action='store_true');return p.parse_args()
if __name__=='__main__':
 a=args()
 try:
  if a.command=='progress':progress(a)
  else:run(a)
 except (InterruptedError,KeyboardInterrupt) as ex:log('SAFE STOP',ex);sys.exit(75)
