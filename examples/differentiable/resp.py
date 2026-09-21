"""自洽响应项有多大? 比较 (a) 全自洽重收敛的有限差分 与 (b) 冻结密度的有限差分。"""
import sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/scratch/work/yuc11/altharm_runs/defumat_p7/defumat_src")
import numpy as np, copy
from pathlib import Path
from defumat import Calculator
from defumat.system.kpoints import KPoints
from defumat.workflows.bands import run_bands
PS = Path("/scratch/work/yuc11/altharm_runs/defumat_p7/pseudo_v2")
KTEST = np.array([[0.13, 0.21, 0.07]]); NB = 10
D = 0.004      # 分数坐标位移

def eigen(scf, calc):
    b = run_bands(calc.system, calc.pseudos, scf.density,
                  kpoints=KPoints.from_cartesian(KTEST, weights=np.ones(1)), nbnd=NB,
                  ns=getattr(scf,'ns',None), becsum=getattr(scf,'becsum',()) or ())
    return np.asarray(b.eigenvalues_by_spin)[0,0]

def write_in(frac2, path):
    txt = open("scf_two.in").read().splitlines()
    out=[]
    for ln in txt:
        if ln.strip().startswith("As 0.48"):
            out.append(f"  As {frac2[0]:.10f} {frac2[1]:.10f} {frac2[2]:.10f}")
        else: out.append(ln)
    open(path,"w").write("\n".join(out)+"\n")

base=np.array([0.48,0.51,0.53])
calc0 = Calculator.from_file("scf_two.in", pseudo_dir=PS); scf0 = calc0.get_scf()
print(f"[基准] SCF {scf0.iterations} 步 收敛={scf0.converged} E={scf0.total_energy:.6f} Ry", flush=True)
e0 = eigen(scf0, calc0)

res={}
for sign,lab in ((+1,"p"),(-1,"m")):
    f = base.copy(); f[0] += sign*D
    write_in(f, f"scf_{lab}.in")
    c = Calculator.from_file(f"scf_{lab}.in", pseudo_dir=PS)
    s = c.get_scf()                                    # (a) 全自洽
    res[("scf",lab)] = eigen(s, c)
    res[("frz",lab)] = eigen(scf0, c)                  # (b) 用基准密度, 只动了原子
    print(f"[{lab}] 全自洽 SCF {s.iterations} 步 收敛={s.converged}", flush=True)

d_scf = (res[("scf","p")]-res[("scf","m")])/(2*D)
d_frz = (res[("frz","p")]-res[("frz","m")])/(2*D)
print(f"\nd(eps)/d(x_2), 分数坐标, Ry:")
print(f"{'带':>3} {'全自洽 (真值)':>15} {'冻结密度':>14} {'自洽响应项':>14} {'占比':>8}")
print("-"*62)
for n in range(8):
    diff=d_scf[n]-d_frz[n]
    frac = abs(diff)/max(abs(d_scf[n]),1e-9)*100
    print(f"{n:3d} {d_scf[n]:15.6f} {d_frz[n]:14.6f} {diff:14.6f} {frac:7.1f}%")
m=[abs(d_scf[n]-d_frz[n])/max(abs(d_scf[n]),1e-9) for n in range(8)]
print(f"\n>>> 自洽响应项占真值的比例: 中位 {np.median(m)*100:.1f}%, 最大 {max(m)*100:.1f}%")
