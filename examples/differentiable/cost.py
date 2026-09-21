import sys, time, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/scratch/work/yuc11/altharm_runs/defumat_p7/defumat_src")
import numpy as np
from pathlib import Path
from defumat import Calculator
from defumat.system.kpoints import KPoints
from defumat.response.velocity import band_velocities
from defumat.workflows.bands import run_bands
calc = Calculator.from_file("scf_tiny.in", pseudo_dir=Path("/scratch/work/yuc11/altharm_runs/defumat_p7/pseudo_v2"))
scf = calc.get_scf()
NB=16; rng=np.random.default_rng(3)
for nk in (8, 24):
    ks = rng.uniform(-0.4,0.4,(nk,3))
    kp = KPoints.from_cartesian(ks, weights=np.ones(nk))
    t=time.time(); run_bands(calc.system, calc.pseudos, scf.density, kpoints=kp, nbnd=NB,
                             ns=getattr(scf,'ns',None), becsum=getattr(scf,'becsum',()) or ()); tb=time.time()-t
    t=time.time(); band_velocities(calc, scf, kpoints=kp, nbnd=NB); tv=time.time()-t
    print(f"nk={nk:3d}  bands {tb:6.2f}s   bands+velocities {tv:6.2f}s   倍数 {tv/tb:.2f}x")
