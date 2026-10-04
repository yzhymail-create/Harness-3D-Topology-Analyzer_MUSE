#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""线束拓扑校验逻辑 v1 (2026-10-04)
对 topology.json 做分级检查:
  ERROR - 硬性规则, 不过则拓扑不可用, 必须修代码
  WARN  - 预警, 需人工复核几何是否合理
用法: python3 check_topology.py <topology.json>
"""
import json
import math
import sys
from collections import Counter


def check(topo):
    errors, warns = [], []
    nodes = topo.get("nodes", [])
    segs = topo.get("segments", [])
    runs = topo.get("runs", [])
    nc = {n["code"]: n for n in nodes}
    sg = {s["seg"]: s for s in segs}

    # ---- A. 结构完整性 (ERROR) ----
    for s in segs:  # A1/A2: 端点存在 + 双向一致
        for e in ("na", "nb"):
            if s[e] not in nc:
                errors.append(f"A1 线段{s['seg']}端点{s[e]}无对应节点")
            elif s["seg"] not in nc[s[e]].get("segments", []):
                errors.append(f"A2 线段{s['seg']}与节点{s[e]}互相引用不一致")
        if not (s.get("length") or 0) > 0:  # A4
            errors.append(f"A4 线段{s['seg']}长度非正")
    for n in nodes:  # A3
        if n.get("degree", 0) < 1:
            errors.append(f"A3 节点{n['code']}孤立(degree=0)")
        xyz = n.get("xyz") or []
        if len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):  # E1
            warns.append(f"E1 节点{n['code']}坐标异常")

    # ---- B. 用户规则 (ERROR) ----
    for n in nodes:
        c, deg = n["code"], n.get("degree", 0)
        if c.startswith("CON") and deg != 1:  # B1
            errors.append(f"B1 连接器{c}度数={deg}(规则: 只能连一个线段端点)")
        if (c.startswith("BN") or n.get("type") == "fork") and deg < 3:  # B2
            errors.append(f"B2 分支点{c}度数={deg}(规则: >=3)")
        if (c.startswith("N") or n.get("type") == "hanging") and deg != 1:  # B3
            warns.append(f"B3 自由端{c}度数={deg}")

    # ---- C. 走线 (ERROR) ----
    for r in runs:  # C1: 连续性
        ss = r.get("segs", [])
        for a, b in zip(ss, ss[1:]):
            if a not in sg or b not in sg:
                errors.append(f"C1 走线{r['run']}含未知线段")
                break
            sa, sb = sg[a], sg[b]
            if not ({sa["na"], sa["nb"]} & {sb["na"], sb["nb"]}):
                errors.append(f"C1 走线{r['run']}中{a}-{b}不共享节点")
                break
    in_run = Counter()
    for r in runs:  # C2: 覆盖
        for s in r.get("segs", []):
            in_run[s] += 1
    for s in segs:
        if in_run[s["seg"]] != 1:
            warns.append(f"C2 线段{s['seg']}属{in_run[s['seg']]}条走线")

    # ---- D. 分类合理性 ----
    for n in nodes:  # D1: 卡扣度数
        c, deg = n["code"], n.get("degree", 0)
        if c.startswith("CLP") and deg == 1:
            # 用户规则: 卡扣不能在顶端, degree=1 的卡扣是分类错误
            errors.append(f"D1 卡扣{c}在顶端degree=1(应为连接器), 分类错误")
        elif c.startswith("CLP") and deg >= 3:
            warns.append(f"D1 卡扣{c}度数={deg}(交汇处, 需人工确认)")
    branches = {s["branch"] for s in segs}  # D3: 分支覆盖
    for b in topo.get("branches", []):
        if b.get("key") not in branches:
            warns.append(f"D3 分支{b.get('key')}无对应线段")

    return errors, warns


def main():
    path = sys.argv[1]
    topo = json.load(open(path, encoding="utf-8"))
    errors, warns = check(topo)
    n_seg, n_node, n_run = len(topo.get("segments", [])), len(topo.get("nodes", [])), len(topo.get("runs", []))
    print(f"拓扑: {n_seg}线段/{n_node}节点/{n_run}走线")
    print(f"ERROR: {len(errors)}  WARN: {len(warns)}")
    for e in errors:
        print("  [ERROR]", e)
    for w in warns:
        print("  [WARN] ", w)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
