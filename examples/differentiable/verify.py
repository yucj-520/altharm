import sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/scratch/work/yuc11/altharm_runs/defumat_p7/defumat_src")
import numpy as np
from pathlib import Path
from defumat import Calculator
from defumat.system.kpoints import KPoints
from defumat.response.velocity import band_velocities
from defumat.workflows.bands import run_bands

BOHR=0.52917720859; ALAT_A=2.9; ALAT=ALAT_A/BOHR; TPIBA=2*np.pi/ALAT
calc = Calculator.from_file("scf_tiny.in", pseudo_dir=Path("/scratch/work/yuc11/altharm_runs/defumat_p7/pseudo_v2"))
scf = calc.get_scf()
print(f"alat={ALAT:.5f} bohr, 2pi/alat={TPIBA:.6f} bohr^-1\n")
NB=8; rng=np.random.default_rng(7)
worst=0.0
for trial in range(3):
    k0 = rng.uniform(-0.3,0.3,3)          # 一般位置, 避开简并
    bv = band_velocities(calc, scf, kpoints=KPoints.from_cartesian(k0[None,:],weights=np.ones(1)), nbnd=NB)
    v = np.asarray(bv.velocities)
    if v.ndim==4: v=v[0]
    for axis,name in ((0,'x'),(1,'y'),(2,'z')):
        d=5e-4; e=np.zeros(3); e[axis]=d
        b = run_bands(calc.system, calc.pseudos, scf.density,
                      kpoints=KPoints.from_cartesian(np.array([k0+e,k0-e]),weights=np.ones(2)),
                      nbnd=NB, ns=getattr(scf,'ns',None), becsum=getattr(scf,'becsum',()) or ())
        ev=np.asarray(b.eigenvalues_by_spin)[0]
        fd=(ev[0]-ev[1])/(2*d)                     # Ry per tpiba
        auto=v[0,:,axis]*TPIBA                     # Ry bohr -> Ry per tpiba
        rel=np.abs(auto-fd)/np.maximum(np.abs(fd),1e-6)
        worst=max(worst,rel[:6].max())
        if trial==0 and axis==0:
            print("  带   autodiff      有限差分      相对差")
            for n in range(6):
                print(f"  {n:2d}  {auto[n]:+.6f}  {fd[n]:+.6f}   {rel[n]:.2e}")
    print(f"trial {trial}: k={np.round(k0,3)}  三个方向前六带最大相对差 = {worst:.2e}")
print(f"\n>>> autodiff 与有限差分最大相对偏差: {worst:.2e}")
