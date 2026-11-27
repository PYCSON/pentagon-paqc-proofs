#!/usr/bin/env python3
"""Full qtable fibre/closed stage (shared checkpoint engine, default 20 workers)."""
import os,sys,pathlib
os.environ["PYTHONOPTIMIZE"]="0"
p=pathlib.Path(__file__).resolve().parent
os.execv(sys.executable,[sys.executable,str(p/'code/driver.py'),'run','--stage','q_other',*sys.argv[1:]])
