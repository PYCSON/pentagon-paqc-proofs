#!/usr/bin/env python3
import sys,pathlib
if sys.flags.optimize: raise SystemExit("Do not run the strict verifier with Python -O")
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent/'code'))
from replay_all_impl import main
if __name__=='__main__':main()
