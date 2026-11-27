#!/usr/bin/env python3
"""Exact witness binding, tamper and resume tests; only independent code."""
import copy,json,pathlib,sys,tempfile,time
from driver import enumerate_stage,identity,canonical,load_shard,guard
from replay_all_impl import verify_one
from types import SimpleNamespace
P=pathlib.Path
if __name__=='__main__':
 data=P(sys.argv[1]);results=P(sys.argv[2]);
 rows=[json.loads(s) for s in (results/'near/000000.jsonl').open()];r=next(x for x in rows if x['status']=='PROVED');t=next(t for t in enumerate_stage(data,'near') if t['id']==r['id'])
 assert verify_one((t,r))[1]=='PROVED';print('PRISTINE PASS',r['id'],flush=True)
 def reject(name,r1,t1=t):
  assert verify_one((t1,r1))[1]=='REJECT',name;print('MUST_REJECT PASS',name,flush=True)
 x=copy.deepcopy(r);x['task_sha256']='0'*64;reject('forged identity',x)
 x=copy.deepcopy(r);del x['witness'];reject('missing witness',x)
 x=copy.deepcopy(r);x['status']='UNPROVED';reject('nonpassing record with witness',x)
 x=copy.deepcopy(r);w=x['witness'];
 while 'lg' not in w:w=w['left']
 w['lg']['digest']='0'*64;reject('LG polynomial digest',x)
 x=copy.deepcopy(r);w=x['witness'];
 while 'lg' not in w:w=w['left']
 w['lg']['L']='0';reject('weakened target',x)
 x=copy.deepcopy(r);w=x['witness'];
 while 'lg' not in w:w=w['left']
 w['lg']['vertices']=[];reject('geometry override',x)
 x=copy.deepcopy(r);w=x['witness'];
 while 'lg' not in w:w=w['left']
 w['lg']['lo'][0]='999';reject('wrong child box',x)
 with tempfile.TemporaryDirectory() as td:
  p=P(td)/'journal';line=(canonical(r)+'\n').encode();p.write_bytes(line+b'{"broken');got=load_shard(p,{t['id']:t},repair=True,metadata_only=True);assert len(got)==1 and p.read_bytes()==line;print('CHECKPOINT truncated tail recovered',flush=True)
  p.write_bytes(line+line)
  try:load_shard(p,{t['id']:t},repair=True);raise AssertionError('duplicate accepted')
  except ValueError:print('MUST_REJECT PASS duplicate record',flush=True)
  try:guard(P(td),SimpleNamespace(disk_min=1e9,mem_min=2));raise AssertionError('disk guard did not stop')
  except InterruptedError as ex: print('DISK GUARD PASS',ex,flush=True)
 print('CONTROLS PASS',flush=True)
