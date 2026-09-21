# 当前状态 — 2026-09-21

会话中断后从这里接上。**仓库和集群上的东西都不受会话关闭影响**；受影响的只有会话里的后台等待和 subagent。

## 结论口径（会变的只有这一节）

**参考值**：QE 24×24×8 + `degauss = 0.0150`，band 44（0-based），152 点圆周，k_z = 0.25·b₃。

| l | 3 | 9 | 15 | 21 | 27 | 33 | 39 | 45 |
|---|---|---|---|---|---|---|---|---|
| meV | 64.293 | 43.447 | 3.603 | 8.662 | 1.876 | 2.583 | 1.182 | 0.985 |

噪声底 23.61 μeV。日常计算用 12×12×4 @ 0.0150（离参考值中位 1.2%、最差 2.1%）。

**旧口径（degauss 0.0037）全部作废**，且不能靠乘系数修正——见 PLAN.md 的 P6b。

## 正在跑的作业（会话关掉也继续跑）

| 作业号 | 名字 | 内容 | 回来后怎么看 |
|---|---|---|---|
| 20367432 | dfm_corr | defumat 在订正设置下重跑，补三方表 | `altharm_runs/defumat_p7/dfm_corr-*.out`，产物 `defumat_corr.pkl` |
| 20367531/5/6/7/8 | g1e3 g1e4 g1e2 gz1e3 gxy1e3 | 几何扰动敏感性（subagent 投的） | `altharm_runs/` 下对应目录 |

`squeue -u $USER` 看队列。

## 会话关闭会丢掉的

- 后台等待（作业结束时的自动通知）——作业本身不受影响，回来手动看就行
- **几何精度那个 subagent**——它投的 5 个 SLURM 作业会继续跑完，但它本身的分析和写测试的进度会中断。`tests/test_geometry_precision.py` 是它未完成的产物，尚未提交。回来后要么让我接手分析那批结果，要么重新派一个。

## 回来后的下一步

1. 看 dfm_corr 和几何扰动的结果，补进 PLAN.md 和 PPT（P6/P7 两页的数字仍是旧口径）
2. **Hubbard U 敏感性从未测过**——U = 3.0 eV 是课题组选的，而噪声底看不见它。这是我建议的下一个优先项
3. 仓库仍未推到 GitHub：`git remote add origin https://github.com/yucj-520/altharm.git && git push -u origin master`

## 磁盘

配额 112 G / 500 G。`altharm_runs` 114 M。清理规则：`AB.wfc*`、`AB.save/wfc*.dat`、`AB.hub*` 可删，**但 harmwave.x 需要波函数**（见 PLAN.md P9）。
