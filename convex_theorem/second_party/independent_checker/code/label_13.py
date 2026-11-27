#!/usr/bin/env python3
"""Recheck addendum target keys, all matching label records (counts known to be 13 only for selected original rows).
The addendum lacks full box IDs for repeated keys; locate every matching record and report actual counts."""
import json,sys,time,collections,pathlib
from label_runner import recheck
from multiprocessing import Pool
P=pathlib.Path('/mnt/data/p24audit/data/far/out')
entries=[('E3','7_3',[7,3,-1,0],0,'T'),('E3','16_5',[16,5,8,0],3,'T'),('E3','21_11',[21,11,10,1],0,'T'),('E2','15_8',[15,8,11,1],3,'TX'),('E2','1_4',[1,4,9,0],0,'T'),('E2','-9_11',[-9,11,7,0],0,'T'),('E3','15_7',[15,7,5,0],3,'T'),('E3','16_6',[16,6,7,1],0,'T'),('E2','8_4',[8,4,1,1],0,'T'),('E2','-4_5',[-4,5,8,1],3,'T')]
def work(r):
 st=time.monotonic();chart,block,line,k,j,lo,hi,lab=r
 try:s,sp=recheck(chart,j,lo,hi,lab,0,1);return dict(chart=chart,block=block,line=line,key=k,j=j,lo=lo,hi=hi,label=lab,status=s or 'UNPROVED',split=sp,seconds=time.monotonic()-st)
 except Exception as e:return dict(chart=chart,block=block,line=line,key=k,j=j,label=lab,status='ERROR',reason=repr(e))
if __name__=='__main__':
 rows=[]
 for c,b,k,j,l in entries:
  for i,line in enumerate(open(P/c/'blocks'/b/'leaves.jsonl'),1):
   kk,jj,lo,hi,ll=json.loads(line)
   if [kk,jj,ll]==[k,j,l]:rows.append((c,b,i,k,j,lo,hi,l))
 print('MATCHING',len(rows),flush=True);ct=collections.Counter()
 with Pool(1) as p,open('logs/B_addendum_keys.jsonl','w',buffering=1) as f:
  for r in p.imap_unordered(work,rows,chunksize=1):ct[r['status'],r.get('split')]+=1;f.write(json.dumps(r)+'\n');print(r['chart'],r['block'],r['line'],r['status'],r.get('split'),flush=True)
 print('SUMMARY',ct,flush=True)
