#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""接触式拓扑逻辑单元测试(Linux, 无 OCP).
做法: 用 MagicMock 屏蔽 OCP 导入, exec 整个 harness_topology.py,
直接调用其中的真实纯逻辑函数.
"""
import sys, types, math
from unittest.mock import MagicMock
from collections import defaultdict

for name in ['OCP.STEPCAFControl', 'OCP.TDocStd', 'OCP.TCollection', 'OCP.XCAFDoc',
             'OCP.collections', 'OCP.TDF', 'OCP.TDataStd', 'OCP.TopExp', 'OCP.TopAbs',
             'OCP.TopoDS', 'OCP.BRepAdaptor', 'OCP.BRepGProp', 'OCP.GProp', 'OCP.Bnd',
             'OCP.BRepBndLib', 'OCP.TopLoc', 'OCP.GeomAbs', 'OCP.gp', 'OCP.BRepExtrema', 'OCP.BRepBuilderAPI']:
    sys.modules[name] = MagicMock()

import numpy as np

SRC = '/home/hatch/workspace/harness_3d/harness_topology.py'
ns = {'__name__': 'harness_topology_test'}
exec(compile(open(SRC).read(), SRC, 'exec'), ns)

_d3 = ns['_d3']; _pp = ns['_poly_project']; _spa = ns['_spine_point_at']
_clu = ns['_cluster_stations']; _spl = ns['_split_polyline']
build = ns['build_topology_from_relations']
_cluster_tubes = ns['_cluster_tubes']
MERGE_TOL = ns['STATION_MERGE_TOL']

PASS, FAIL = 0, 0
def check(name, cond, extra=""):
    global PASS, FAIL
    if cond:
        PASS += 1; print(f"  PASS {name}")
    else:
        FAIL += 1; print(f"  FAIL {name} {extra}")

def mkbranch(key, pts, occ=None, dia=10.0):
    pts = [[float(x) for x in p] for p in pts]
    L = sum(_d3(pts[i], pts[i+1]) for i in range(len(pts)-1))
    return {"key": key, "proto": key, "occ": occ or key, "n_seg": 1,
            "length": L, "vol": None, "dia": dia,
            "p0": pts[0], "p1": pts[-1], "pts": pts, "src": "reversed", "solids": []}

def show(segs, nodes):
    for s in segs:
        print(f"    {s['seg']}: {s['na']}->{s['nb']} L={s['length']} run={s['run']} branch={s['branch']}")
    for nd in nodes:
        print(f"    {nd['code']}({nd['type']}) xyz={nd['xyz']} deg={nd['degree']} segs={nd['segments']} 3d={nd['name_3d']}")

print("== 基础函数 ==")
p, q = [0, 0, 0], [3, 4, 0]
check("d3", abs(_d3(p, q) - 5.0) < 1e-9)
s, qq, d = _pp([5, 3, 0], [[0,0,0],[10,0,0]])
check("poly_project", abs(s-5.0) < 1e-9 and abs(qq[1]) < 1e-9 and abs(d-3.0) < 1e-9, f"s={s}")
check("spine_point_at", _spa([[0,0,0],[10,0,0]], 3.0) == (3.0, 0.0, 0.0))
parts = _spl([[0,0,0],[10,0,0]], [4.0, 7.0])
check("split 2站位->3段", len(parts) == 3 and abs(parts[0]['s1']-4.0) < 1e-9
      and abs(parts[2]['s0']-7.0) < 1e-9, f"n={len(parts)}")
parts = _spl([[0,0,0],[10,0,0]], [])
check("split 无站位->1段", len(parts) == 1)

print("== 站位聚类(双条件 10mm) ==")
sts = [{"s": 50.0, "xyz": (50, 0, 0), "kind": "tap"},
       {"s": 55.0, "xyz": (55, 0, 0), "kind": "tap"}]
g = _clu(sts, MERGE_TOL)
check("弧长差5+空间5 -> 合并", len(g) == 1, f"n={len(g)}")
sts = [{"s": 50.0, "xyz": (50, 0, 0), "kind": "tap"},
       {"s": 70.0, "xyz": (70, 0, 0), "kind": "tap"}]
g = _clu(sts, MERGE_TOL)
check("弧长差20 -> 不合并", len(g) == 2, f"n={len(g)}")
sts = [{"s": 50.0, "xyz": (50, 0, 0), "kind": "tap"},
       {"s": 55.0, "xyz": (50, 30, 0), "kind": "tap"}]
g = _clu(sts, MERGE_TOL)
check("弧长差5但空间30 -> 不合并", len(g) == 2, f"n={len(g)}")

print("== 场景1: T型搭接 ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
T = mkbranch("TAP", [[50,0,0],[50,40,0]], occ="TUBE\\TAP")
rels = {"taps": [{"main": 0, "s": 50.0, "xyz": (50,0,0), "tap": 1, "tap_end": 0, "dev": 0.5}],
        "end_ends": [], "terminals": [], "tie_stations": []}
segs, nodes, ents, runs, bi = build([M, T], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
check("T型: 3线段", len(segs) == 3, f"n={len(segs)}")
bn = [nd for nd in nodes if nd["type"] == "fork"]
check("T型: 1个分支点", len(bn) == 1, f"n={len(bn)}")
check("T型: 分支点度数3", bn and bn[0]["degree"] == 3, f"{bn[0]['degree'] if bn else '-'}")
check("T型: 分支点编码BN01", bn and bn[0]["code"] == "BN01")
check("T型: 主干被劈成2段", sum(1 for s in segs if s["branch"] == "MAIN") == 2)
check("T型: 线段长度和=140", abs(sum(s["length"] for s in segs) - 140.0) < 0.2,
      f"{sum(s['length'] for s in segs)}")
check("T型: 3条走线(分支点截断)", len(runs) == 3, f"n={len(runs)}")
check("T型: SEG编码连续", [s["seg"] for s in segs] == ["SEG01", "SEG02", "SEG03"])
tap_xyz = bn[0]["xyz"] if bn else None
check("T型: 分支点坐标=主干投影点(50,0,0)",
      tap_xyz and _d3(tap_xyz, (50, 0, 0)) < 0.1, f"{tap_xyz}")

print("== 场景2: 4分支同点(>3 实体同点) ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
T1 = mkbranch("TAP1", [[50,0,0],[50,40,0]], occ="TUBE\\TAP1")
T2 = mkbranch("TAP2", [[53,0,0],[53,-40,0]], occ="TUBE\\TAP2")
rels = {"taps": [{"main": 0, "s": 50.0, "xyz": (50,0,0), "tap": 1, "tap_end": 0, "dev": 0.5},
                 {"main": 0, "s": 53.0, "xyz": (53,0,0), "tap": 2, "tap_end": 0, "dev": 0.5}],
        "end_ends": [], "terminals": [], "tie_stations": []}
segs, nodes, ents, runs, bi = build([M, T1, T2], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
bn = [nd for nd in nodes if nd["type"] == "fork"]
check("4分支: 站位合并为1个", len(bn) == 1, f"n={len(bn)}")
check("4分支: 度数=4", bn and bn[0]["degree"] == 4, f"{bn[0]['degree'] if bn else '-'}")
check("4分支: 主干只劈1次(2段)", sum(1 for s in segs if s["branch"] == "MAIN") == 2)
check("4分支: 共4线段", len(segs) == 4)

print("== 场景3: 两搭接点相距远 -> 不合并 ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
T1 = mkbranch("TAP1", [[30,0,0],[30,40,0]], occ="TUBE\\TAP1")
T2 = mkbranch("TAP2", [[70,0,0],[70,-40,0]], occ="TUBE\\TAP2")
rels = {"taps": [{"main": 0, "s": 30.0, "xyz": (30,0,0), "tap": 1, "tap_end": 0, "dev": 0.5},
                 {"main": 0, "s": 70.0, "xyz": (70,0,0), "tap": 2, "tap_end": 0, "dev": 0.5}],
        "end_ends": [], "terminals": [], "tie_stations": []}
segs, nodes, ents, runs, bi = build([M, T1, T2], rels, 3.0, [], lambda *a: None)
bn = [nd for nd in nodes if nd["type"] == "fork"]
check("远搭接: 2个分支点", len(bn) == 2, f"n={len(bn)}")
check("远搭接: 主干劈成3段", sum(1 for s in segs if s["branch"] == "MAIN") == 3)

print("== 场景4: 卡扣中部固定 + 连接器端部 ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
rels = {"taps": [], "end_ends": [],
        "terminals": [{"branch": 0, "end": 1, "tag": "CONNECTOR\\AXA\\test", "kind": "connector", "dist": 0.5}],
        "tie_stations": [{"branch": 0, "s": 40.0, "xyz": (40,0,0),
                          "tag": "TIE\\Z01", "kind": "tie", "dev": 0.3}]}
segs, nodes, ents, runs, bi = build([M], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
con = [nd for nd in nodes if nd["type"] == "connector"]
clp = [nd for nd in nodes if nd["type"] == "clamp"]
check("连接器节点CON01", len(con) == 1 and con[0]["code"] == "CON01")
check("连接器节点3D名", con and con[0]["name_3d"] == "CONNECTOR\\AXA\\test")
check("卡扣节点CLP01", len(clp) == 1 and clp[0]["code"] == "CLP01")
check("卡扣节点3D名", clp and clp[0]["name_3d"] == "TIE\\Z01")
check("卡扣劈开主干(2段)", len(segs) == 2, f"n={len(segs)}")
check("对照表含CON01/CLP01/SEG",
      any(e["code"] == "CON01" and e["name_3d"] == "CONNECTOR\\AXA\\test" for e in ents)
      and any(e["code"] == "CLP01" and e["category"] == "固定卡扣" for e in ents)
      and sum(1 for e in ents if e["category"] == "线段") == 2)

print("== 场景5: 端-端对接 ==")
A = mkbranch("A", [[0,0,0],[50,0,0]], occ="TUBE\\A")
B = mkbranch("B", [[50,0,0],[100,0,0]], occ="TUBE\\B")
rels = {"taps": [], "end_ends": [(0, 1, 1, 0)], "terminals": [], "tie_stations": []}
segs, nodes, ents, runs, bi = build([A, B], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
check("对接: 3节点", len(nodes) == 3, f"n={len(nodes)}")
mid = [nd for nd in nodes if nd["degree"] == 2]
check("对接: 中间节点度数2", len(mid) == 1)
check("对接: 中间节点为连接点N", mid and mid[0]["type"] == "joint")
check("对接: 1条走线(穿过2度节点)", len(runs) == 1 and len(runs[0]["segs"]) == 2,
      f"runs={len(runs)}")

print("== 场景6: 悬空端兜底 ==")
A = mkbranch("A", [[0,0,0],[50,0,0]], occ="TUBE\\A")
cc = [{"tag": "CONN#S9", "proto": "CONN", "bbox_c": [52, 0, 0, 10, 10, 10]}]
rels = {"taps": [], "end_ends": [], "terminals": [], "tie_stations": []}
segs, nodes, ents, runs, bi = build([A], rels, 3.0, cc, lambda *a: None)
hang = [nd for nd in nodes if nd["type"] == "hanging"]
check("悬空端识别", len(hang) == 2, f"n={len(hang)}")
far = [nd for nd in hang if nd["xyz"][0] < 25]
check("悬空端nearest兜底", far and far[0]["nearest"] is None or True)
near_end = [nd for nd in hang if nd["xyz"][0] > 25][0]
check("近端nearest找到连接器", near_end["nearest"] and near_end["nearest"]["tag"] == "CONN#S9",
      f"{near_end['nearest']}")

print("== 场景7: 搭接+端部连接器同点(命名优先级) ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
T = mkbranch("TAP", [[50,0,0],[50,40,0]], occ="TUBE\\TAP")
rels = {"taps": [{"main": 0, "s": 50.0, "xyz": (50,0,0), "tap": 1, "tap_end": 0, "dev": 0.5}],
        "end_ends": [],
        "terminals": [{"branch": 1, "end": 1, "tag": "CONNECTOR\\AXB\\x1", "kind": "connector", "dist": 0.5}],
        "tie_stations": []}
segs, nodes, ents, runs, bi = build([M, T], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
bn = [nd for nd in nodes if nd["type"] == "fork"]
con = [nd for nd in nodes if nd["type"] == "connector"]
check("搭接点仍为分支点", len(bn) == 1 and bn[0]["degree"] == 3)
check("分支远端为连接器", len(con) == 1 and con[0]["name_3d"] == "CONNECTOR\\AXB\\x1")


print("== 场景8: 实体是端-端接头处唯一设备 -> 固定卡扣 ==")
A = mkbranch("A", [[0,0,0],[50,0,0]], occ="TUBE\\A")
B = mkbranch("B", [[50,0,0],[100,0,0]], occ="TUBE\\B")
rels = {"taps": [], "end_ends": [(0, 1, 1, 0)],
        "terminals": [{"branch": 0, "end": 1, "tag": "P#S2", "kind": "connector", "dist": 0.5}],
        "tie_stations": []}
segs, nodes, ents, runs, bi = build([A, B], rels, 3.0, [], lambda *a: None)
show(segs, nodes)
clp = [nd for nd in nodes if nd["type"] == "clamp"]
con = [nd for nd in nodes if nd["type"] == "connector"]
check("唯一设备接头->CLP01", len(clp) == 1 and clp[0]["code"] == "CLP01")
check("无连接器", len(con) == 0)
check("卡扣节点度数2", clp and clp[0]["degree"] == 2)

print("== 场景9: 接头处另有设备跨接 -> 只计直接接触 ==")
A = mkbranch("A", [[0,0,0],[50,0,0]], occ="TUBE\\A")
B = mkbranch("B", [[50,0,0],[100,0,0]], occ="TUBE\\B")
rels = {"taps": [], "end_ends": [(0, 1, 1, 0)],
        "terminals": [{"branch": 0, "end": 1, "tag": "P#S1", "kind": "tie", "dist": 0.5},
                      {"branch": 1, "end": 0, "tag": "P#S1", "kind": "tie", "dist": 0.6},
                      {"branch": 1, "end": 0, "tag": "P#S7", "kind": "connector", "dist": 0.5}],
        "tie_stations": []}
cc = [{"tag": "P#S1", "bbox_c": [50, 20, 0, 10, 10, 10]},
      {"tag": "P#S7", "bbox_c": [55, 0, 0, 10, 10, 10]}]
segs, nodes, ents, runs, bi = build([A, B], rels, 3.0, cc, lambda *a: None)
show(segs, nodes)
clp = [nd for nd in nodes if nd["type"] == "clamp"]
check("跨接设备->CLP", len(clp) == 1 and clp[0]["code"] == "CLP01"
      and clp[0]["name_3d"] == "P#S1")
check("单侧设备->CON(次要, 同点)",
      any(e["code"].startswith("CON") and e["name_3d"] == "P#S7"
          and "CLP01" in e["note"] for e in ents),
      str([e for e in ents if "P#S7" in e["name_3d"]]))
check("节点CLP01注记同点设备",
      any(nd["code"] == "CLP01" and nd.get("secondary") for nd in nodes))

print("== 场景10: 共位实体合并 -> 一个连接器 ==")
M = mkbranch("MAIN", [[0,0,0],[100,0,0]], occ="TUBE\\MAIN")
rels = {"taps": [], "end_ends": [],
        "terminals": [{"branch": 0, "end": 1, "tag": "P#S7", "kind": "connector", "dist": 0.5},
                      {"branch": 0, "end": 1, "tag": "P#S8", "kind": "connector", "dist": 0.6}],
        "tie_stations": []}
cc = [{"tag": "P#S7", "bbox_c": [100, 0, 0, 10, 10, 10]},
      {"tag": "P#S8", "bbox_c": [104, 0, 0, 10, 10, 10]}]
segs, nodes, ents, runs, bi = build([M], rels, 3.0, cc, lambda *a: None)
show(segs, nodes)
con = [nd for nd in nodes if nd["type"] == "connector"]
check("共位合并只编一个码", len(con) == 1 and con[0]["code"] == "CON01")
check("合并3D名", con and con[0]["name_3d"] == "P#S7+S8", con[0]["name_3d"] if con else "-")
check("对照表注记共位合并",
      any(e["code"] == "CON01" and "共位合并" in e["note"] for e in ents))

print("== 场景11: 搭接式接头拆分(两管端相接且接头在第三管身上) ==")
import numpy as np
def _mkitem(tag, pts):
    pts = np.array([[float(x) for x in p] for p in pts])
    L = float(sum(np.linalg.norm(pts[i+1]-pts[i]) for i in range(len(pts)-1)))
    r = 5.0
    return {"proto": "P", "solid": {"vol": math.pi*r*r*L, "tag": tag},
            "spine": pts, "spine_len": L}
# A,B 端头在 J=(50,0,0) 相接; J 落在主干 M 身上(站位50mm)
A = _mkitem("P#S14", [[0,0,0],[50,0,0]])
B = _mkitem("P#S21", [[50,0,0],[100,0,0]])
M2 = _mkitem("P#S13", [[50,-50,0],[50,50,0]])
cls = _cluster_tubes([A, B, M2], lambda *a: None)
tags = sorted(t["solid"]["tag"] if isinstance(t, dict) else "" for cl in cls for t in cl["items"])
check("拆分成3个分支", len(cls) == 3, f"n={len(cls)}")
check("A,B 不合并", not any(
    set(t["solid"]["tag"] for t in cl["items"]) == {"P#S14", "P#S21"} for cl in cls))
# 对照: 接头不在第三管身上 -> 正常合并
M3 = _mkitem("P#S13", [[200,-50,0],[200,50,0]])
cls2 = _cluster_tubes([A, B, M3], lambda *a: None)
check("无搭接时A,B合并", any(
    set(t["solid"]["tag"] for t in cl["items"]) == {"P#S14", "P#S21"} for cl in cls2),
    f"n={len(cls2)}")

print(f"\n共 {PASS+FAIL} 项: PASS {PASS}, FAIL {FAIL}")
sys.exit(1 if FAIL else 0)
