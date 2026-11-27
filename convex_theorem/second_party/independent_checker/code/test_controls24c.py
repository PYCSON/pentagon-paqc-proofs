#!/usr/bin/env python3
"""Private test on authentic frozen sample and controlled corruptions; never mutates source."""
import argparse,copy,json,tempfile,fcntl,subprocess,sys
from pathlib import Path
import driver as D
import replay24c as R

def run():
 p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('data_zip',type=Path);a=p.parse_args();R.check_code();plan=R.load_plan(a.results)
 with tempfile.TemporaryDirectory() as td:
  root=Path(td);data=D.prepare(a.data_zip,root);rows=[json.loads(x) for x in (a.results/'q_lg/000000.jsonl').open()];r=next(v for v in rows if v['status']=='PROVED');t=next(x for x in D.task_iter(data,'q_lg',plan['config']['sample']) if x['id']==r['id']); assert R.verify_one((t,r,False,64))['status']=='PROVED';print('PRISTINE PASS',r['id'],flush=True)
  def leaf(w):
   while 'lg' not in w:w=w['left']
   return w['lg']
  cases=[('identity',lambda x:x.update(task_sha256='0'*64)),('missing_witness',lambda x:x.pop('witness')),('nonpassing_with_witness',lambda x:x.update(status='UNPROVED')),('digest',lambda x:leaf(x['witness']).update(digest='0'*64)),('target',lambda x:leaf(x['witness']).update(L='0')),('geometry_override',lambda x:leaf(x['witness']).update(vertices=[])),('child_box',lambda x:leaf(x['witness'])['lo'].__setitem__(0,'999')),('N0',lambda x:leaf(x['witness']).update(N=0)),('N13',lambda x:leaf(x['witness']).update(N=13)),('N_mismatch_main',lambda x:leaf(x['witness']).update(N=1 if t['N']!=1 else 2)),('family',lambda x:leaf(x['witness']).update(u=[[],[]],chi=[[],[]])),('dimension',lambda x:leaf(x['witness'])['lo'].append('0'))]
  for name,change in cases:
   x=copy.deepcopy(r);change(x);z=R.verify_one((t,x,False,64));assert z['status']=='REJECT',(name,z);print('MUST_REJECT PASS',name,z['reason'],flush=True)
  from types import SimpleNamespace
  path=root/'test.jsonl';path.write_text(D.canonical(r)+'\n'+D.canonical(r)+'\n')
  try:D.load_shard(path,{t['id']:t});raise AssertionError('duplicate accepted')
  except ValueError:print('MUST_REJECT PASS duplicate',flush=True)
  j=root/'journal';R.append(j,dict(shard=1,n=1,binding={'sources':{'main':'x'}}));b=j.read_bytes()
  with j.open('ab') as f:f.write(b'{"partial')
  got=R.journal(j,True,metadata=True);assert j.read_bytes()==b and got[1]['n']==1;print('JOURNAL TAIL RECOVERY PASS',flush=True)
  try:D.guard(root,SimpleNamespace(disk_min=1e9,mem_min=0));raise AssertionError('disk guard')
  except InterruptedError:print('DISK GUARD PASS',flush=True)
  try:D.guard(root,SimpleNamespace(disk_min=0,mem_min=1e9));raise AssertionError('memory guard')
  except InterruptedError:print('MEMORY GUARD PASS',flush=True)
  (root/'STOP').touch()
  try:D.guard(root,SimpleNamespace(disk_min=0,mem_min=0));raise AssertionError('STOP guard')
  except InterruptedError:print('STOP GUARD PASS',flush=True)
  # Hold the actual shared lock and run the CLI; it must stop before any replay/source processing.
  with (root/'run.lock').open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
   z=subprocess.run([sys.executable,'-S',str(Path(R.__file__)),'--work',str(root),'--stage','summary','--disk-min','0','--mem-min','0'],capture_output=True,text=True)
   assert z.returncode==75,z;print('SHARED LOCK PASS',z.stdout.strip(),flush=True)
 print('CONTROLS PASS',flush=True)
if __name__=='__main__':run()
