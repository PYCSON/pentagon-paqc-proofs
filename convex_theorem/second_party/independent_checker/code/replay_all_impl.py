#!/usr/bin/env python3
"""Trusted-side independent verifier. Never imports anything from submitted results.
Reconstructs original task identity and exact coverage from the SHA-locked frozen input.
"""
import argparse,collections,itertools,json,multiprocessing as mp,os,pathlib,re,subprocess,sys,tempfile,time
from fractions import Fraction as Q
from driver import (BASE,STAGES,LG,ZIP_SHA,VERSION,CHUNK,MAX_ROW,code_id,sha,identity,prepare,task_iter,enumerate_stage,chunks,limits,load_shard)
from proofs24b import lg_bound,discrete_bound
P=pathlib.Path

def verify_manifest(root):
 m=root/'SHA256SUMS';seen=set()
 for line in m.read_text().splitlines():
  if not re.fullmatch(r'[0-9a-f]{64}  .+',line):raise ValueError('bad manifest syntax')
  digest,rel=line.split('  ',1);p=P(rel)
  if p.is_absolute() or '..' in p.parts or rel in seen or rel=='SHA256SUMS':raise ValueError('unsafe or duplicate manifest path')
  f=root/p
  if not f.is_file() or any(x.is_symlink() for x in (f,*[x for x in f.parents if x!=root and root in x.parents])) or sha(f)!=digest:raise ValueError('manifest mismatch '+rel)
  seen.add(rel)
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
 if actual!=seen|{'SHA256SUMS'}:raise ValueError('manifest paths differ '+str(actual^(seen|{'SHA256SUMS'})))
 return seen

def verify_one(args):
 t,r=args
 if r['task_sha256']!=identity(t) or r['id']!=t['id']:return t['id'],'REJECT','frozen identity'
 s=r['status'];w=r.get('witness')
 if s in ('UNPROVED','ERROR'):
  return t['id'],s if w is None else 'REJECT','nonpassing record'
 if w is None:return t['id'],'REJECT','missing witness'
 try:
  if t['stage'] in LG:actual='PROVED' if lg_bound(t,w) else 'REJECT'
  else:actual=discrete_bound(t,w,Q(t.get('table','0')))
  return t['id'],actual if actual==s else 'REJECT','' if actual==s else 'witness/status mismatch '+str(actual)
 except Exception as e:return t['id'],'REJECT',repr(e)[:400]

def geometry(data,work,jobs):
 # Rebuild exact structure and chart root domain, including near robust centre and qtable coverage.
 scripts=BASE/'code';sp=work/'structure.jsonl';qp=work/'q_analytic.jsonl';logs={}
 for name,args in [('structure',[str(data/'out'),str(data),'--jobs',str(jobs),'--out',str(sp)]),('near_robust',[str(data)]),('q_audit',[str(data/'out'),'--out',str(qp)])]:
  r=subprocess.run([sys.executable,'-S',str(scripts/(name+'.py')),*args],capture_output=True,text=True)
  print('GEOMETRY',name,'exit',r.returncode,flush=True);print(r.stdout,flush=True)
  if r.stderr:print(r.stderr,flush=True)
  if r.returncode:raise ValueError('geometry checker failure '+name)
  logs[name]=r.stdout
 rows=[json.loads(s) for s in sp.open()]
 if any(r.get('errors') or r.get('missing') for r in rows) or len([r for r in rows if 'block' in r])!=1260 or len([r for r in rows if 'needed' in r])!=3 or len([r for r in rows if 'near' in r])!=1:raise ValueError('structure coverage mismatch')
 if 'ROOTS MATCH True 57093 BAD 0' not in logs['q_audit'] or 'BAD_LIST []' not in logs['q_audit']:raise ValueError('q table geometry mismatch')
 # Exact LB binding, not just equality of key sets.
 roots=json.loads((data/'out/std/roots.json').read_text())['blocks'];d={}
 for rs in roots.values():
  for k,lo,hi,L in rs:
   k=tuple(k)
   if k in d:raise ValueError('duplicate std root')
   d[k]=Q(L)
 seen=set()
 for line in (data/'out/qtable/cells.jsonl').open():
  r=json.loads(line);k=tuple(r['key'])
  if k in seen or k not in d or Q(r['LB'])!=d[k]:raise ValueError('table LB binding')
  seen.add(k)
 if seen!=set(d):raise ValueError('missing table roots')
 print('GEOMETRY PASS: root LB links',len(seen),flush=True)
 return True

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=P);p.add_argument('--data-zip',type=P,default=BASE/'PENT24_bundle_b8ab7c1.zip');p.add_argument('--jobs',type=int,default=4);p.add_argument('--skip-structure',action='store_true',help='SAMPLE ONLY; no mathematical sign-off');a=p.parse_args();assert 1<=a.jobs<=10, 'replay max 10 jobs at 1024MiB each'
 root=a.results.resolve();paths=verify_manifest(root);plan=json.loads((root/'plan.json').read_text());cfg=plan['config'];sample=cfg['sample']
 if cfg['version']!=VERSION or cfg['zip_sha256']!=ZIP_SHA or cfg['code']!=code_id() or cfg['chunk']!=CHUNK or type(sample)!=int or sample<0 or cfg['engine'] not in ('scipy','stdlib') or not 1<=cfg['budget']<=64:raise ValueError('unrecognized plan / code identity')
 if not sample and a.skip_structure:raise ValueError('cannot skip geometry for full replay')
 print('REPLAY frozen',ZIP_SHA,'code',cfg['code'],'sample',sample,flush=True);st=time.monotonic();stats={};used={'plan.json','summary.json'}
 with tempfile.TemporaryDirectory(prefix='pent24b-replay-') as td:
  work=P(td);data=prepare(a.data_zip,work);G=False if a.skip_structure else geometry(data,work,a.jobs)
  with mp.get_context('spawn').Pool(a.jobs,initializer=limits,initargs=(1024*1024**2,),maxtasksperchild=128) as pool:
   for stage in STAGES:
    ct=collections.Counter();expected_count=0;full=sum(1 for _ in enumerate_stage(data,stage));assert full==plan['full_totals'][stage]
    for idx,batch in enumerate(chunks(task_iter(data,stage,sample))):
     expected_count+=len(batch);expected={t['id']:t for t in batch};rel=f'{stage}/{idx:06d}.jsonl';path=root/rel
     if not path.exists():ct['MISSING']+=len(batch);continue
     used.add(rel);rr=load_shard(path,expected,False,metadata_only=True);ct['MISSING']+=len(batch)-len(rr)
     # Stream in bounded batches; 128 large JSON witnesses need not all live in RAM.
     with path.open('rb') as f:
      def records():
       for line in f:
        r=json.loads(line);yield expected[r['id']],r
      stream=records()
      while items:=list(itertools.islice(stream,max(1,2*a.jobs))):
       for i,s,msg in pool.imap_unordered(verify_one,items,chunksize=1):
        ct[s]+=1
        if s=='REJECT':print('REJECT',i,msg,flush=True)
    if expected_count!=plan['totals'][stage] or expected_count!=sum(ct.values()):raise ValueError('coverage counter mismatch '+stage)
    stats[stage]=dict(ct);print('REPLAY STAGE',stage,expected_count,dict(ct),flush=True)
 if paths!=used:raise ValueError('unexpected/missing result files: '+str(paths^used))
 # Do not infer validity from summary; use exact regenerated task results. Conditional labels require q full proof and geometry.
 def all_proved(s):return stats[s].get('PROVED',0)==plan['totals'][s]
 q_ok=not sample and G and all(all_proved(s) for s in ('q_lg','q_other'))
 n_ok=not sample and G and all_proved('near')
 f_ok=not sample and G and q_ok and all(all_proved(s) for s in ('far_std','far_E2','far_E3')) and all(stats[s].get('PROVED',0)+stats[s].get('CONDITIONAL_TABLE',0)==plan['totals'][s] for s in ('B_std','B_E2','B_E3'))
 print('SUMMARY',json.dumps(stats,sort_keys=True),'seconds',time.monotonic()-st,flush=True)
 print('q_table_closed',q_ok,'N computational gate',n_ok,'F computational gate',f_ok,'SAMPLE_ONLY',bool(sample),flush=True)
 print('No automatic theorem signature: reviewer must also approve analytic hand-over and scope.',flush=True)
 if any(v.get('REJECT',0) for v in stats.values()):sys.exit(2)
 if not sample and not (n_ok and f_ok):sys.exit(3)
if __name__=='__main__':main()
