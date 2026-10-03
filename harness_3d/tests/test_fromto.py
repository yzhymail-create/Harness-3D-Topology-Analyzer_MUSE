#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fromto_report 纯逻辑单元测试 (无 OCP).
覆盖: 线段聚类(合并/卡扣分段/分支分段/闭环)、端面圆心坐标、
      名称解析、Dijkstra 寻径、PathMarker、断路诊断、原始表读取.
"""
import json
import os
import sys

sys.path.insert(0, '/home/hatch/workspace/harness_3d')
import fromto_report as fr

PASS, FAIL = 0, 0


def check(name, cond, extra=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("  PASS %s" % name)
    else:
        FAIL += 1
        print("  FAIL %s %s" % (name, extra))


def mknode(code, ntype, degree):
    return {"code": code, "type": ntype, "degree": degree,
            "xyz": [0.0, 0.0, 0.0], "segments": [], "entity": None,
            "secondary": []}


def mkseg(seg, na, nb, length, p0, p1, branch="B"):
    return {"seg": seg, "na": na, "nb": nb, "length": float(length),
            "dia": 10.0, "p0": [float(x) for x in p0],
            "p1": [float(x) for x in p1], "branch": branch}


def wire(nodes, segs):
    for n in nodes:
        n["segments"] = [s["seg"] for s in segs
                         if s["na"] == n["code"] or s["nb"] == n["code"]]
    return nodes, segs


# ---------- 1. 普通 joint 两侧合并 ----------
nodes, segs = wire(
    [mknode("CA", "connector", 1), mknode("N1", "joint", 2),
     mknode("CB", "connector", 1)],
    [mkseg("S1", "CA", "N1", 100, [0, 0, 0], [100, 0, 0], "Body.46"),
     mkseg("S2", "N1", "CB", 60, [100, 0, 0], [160, 0, 0], "Body.47")])
curves = fr.cluster_curves(nodes, segs, lambda b: b)
check("joint两侧合并为1段", len(curves) == 1, "got %d" % len(curves))
check("合并长度=160", curves[0]["length"] == 160.0, str(curves[0]["length"]))
check("合并名称", curves[0]["name"] == "Body.46 + Body.47", curves[0]["name"])
check("合并端点", curves[0]["p0"] == [0.0, 0.0, 0.0]
      and curves[0]["p1"] == [160.0, 0.0, 0.0], str(curves[0]["p0"]))

# ---------- 2. 卡扣处保留分段点 ----------
nodes, segs = wire(
    [mknode("CA", "connector", 1), mknode("N1", "joint", 2),
     mknode("CL1", "clamp", 2), mknode("CB", "connector", 1)],
    [mkseg("S1", "CA", "N1", 100, [0, 0, 0], [100, 0, 0], "Body.46"),
     mkseg("S2", "N1", "CL1", 50, [100, 0, 0], [150, 0, 0], "Body.46"),
     mkseg("S3", "CL1", "CB", 80, [150, 0, 0], [230, 0, 0], "Body.47")])
curves = fr.cluster_curves(nodes, segs, lambda b: b)
check("卡扣分段->2段", len(curves) == 2, "got %d" % len(curves))
by_ab = {(c["node_a"], c["node_b"]): c for c in curves}
check("第一段 CA-CL1", ("CA", "CL1") in by_ab or ("CL1", "CA") in by_ab)
c1 = by_ab.get(("CA", "CL1")) or by_ab.get(("CL1", "CA"))
check("第一段长度150", c1["length"] == 150.0, str(c1["length"]))

# ---------- 3. 分支点分段 ----------
nodes, segs = wire(
    [mknode("CA", "connector", 1), mknode("BN1", "fork", 3),
     mknode("CB", "connector", 1), mknode("CC", "connector", 1)],
    [mkseg("S1", "CA", "BN1", 100, [0, 0, 0], [100, 0, 0], "B1"),
     mkseg("S2", "BN1", "CB", 60, [100, 0, 0], [160, 0, 0], "B2"),
     mkseg("S3", "BN1", "CC", 70, [100, 0, 0], [100, 70, 0], "B3")])
curves = fr.cluster_curves(nodes, segs, lambda b: b)
check("分支点分段->3段", len(curves) == 3, "got %d" % len(curves))

# ---------- 4. 全joint闭环不死循环 ----------
nodes, segs = wire(
    [mknode("N1", "joint", 2), mknode("N2", "joint", 2), mknode("N3", "joint", 2)],
    [mkseg("S1", "N1", "N2", 10, [0, 0, 0], [10, 0, 0]),
     mkseg("S2", "N2", "N3", 10, [10, 0, 0], [10, 10, 0]),
     mkseg("S3", "N3", "N1", 10, [10, 10, 0], [0, 0, 0])])
curves = fr.cluster_curves(nodes, segs, lambda b: b)
check("闭环收敛为1段", len(curves) == 1, "got %d" % len(curves))
check("闭环长度30", curves[0]["length"] == 30.0)

# ---------- 5. 端面圆心坐标 ----------
nodes, segs = wire(
    [mknode("CA", "connector", 1), mknode("CL1", "clamp", 2),
     mknode("CB", "connector", 1)],
    [mkseg("S1", "CA", "CL1", 100, [1, 2, 3], [101, 2, 3], "B1"),
     mkseg("S2", "CL1", "CB", 80, [101.5, 2, 3], [181, 2, 3], "B2")])
seg_by_id = {s["seg"]: s for s in segs}
nodes[0]["xyz"] = [50, 50, 50]   # 连接器包围盒中心(故意偏离)
nodes[1]["xyz"] = [101.2, 2.1, 3.0]
check("连接器=线段端面圆心",
      fr.node_face_center(nodes[0], seg_by_id) == [1.0, 2.0, 3.0],
      str(fr.node_face_center(nodes[0], seg_by_id)))
check("卡扣=最近线段端面圆心",
      fr.node_face_center(nodes[1], seg_by_id) == [101.0, 2.0, 3.0],
      str(fr.node_face_center(nodes[1], seg_by_id)))

# ---------- 6. 名称 ----------
nm = {"S7": "CONNECTOR\\BMS.1", "S8": "CONNECTOR\\BMS.2"}
check("多S标签解析", fr.tag_real_names("X#S7+S8", nm)
      == ["CONNECTOR\\BMS.1", "CONNECTOR\\BMS.2"])
check("名称归一化", fr.norm_name("FLOOR_AllCATPart::INLINE-FLH_")
      == "FLOORALLCATPART::INLINE-FLH")

# ---------- 7. 寻径 + PathMarker(含卡扣) ----------
nodes, segs = wire(
    [mknode("CA", "connector", 1), mknode("CL1", "clamp", 2),
     mknode("CB", "connector", 1), mknode("CC", "connector", 1)],
    [mkseg("S1", "CA", "CL1", 150, [0, 0, 0], [150, 0, 0], "Body.46"),
     mkseg("S2", "CL1", "CB", 80, [150, 0, 0], [230, 0, 0], "Body.47"),
     mkseg("S3", "CL1", "CC", 90, [150, 0, 0], [150, 90, 0], "Body.48")])
for n in nodes:
    n["_fx"] = fr.node_face_center(n, {s["seg"]: s for s in segs})
by_code = {n["code"]: n for n in nodes}
name_map = {}
curves = fr.cluster_curves(nodes, segs, lambda b: b)
adj = fr.build_graph(curves)
r = fr.shortest_path(adj, curves, "CA", "CB")
check("CA->CB连通", r is not None)
node_path, curve_path, length = r
check("最短长度230", length == 230.0, str(length))
ca = {"partname": "中控大屏", "name": "中控大屏"}
cb = {"partname": "组合仪表显示屏", "name": "组合仪表显示屏"}
nodes[1]["entity"] = "X#S0"
nm2 = {"S0": "C14\\Tyton-157-00181"}
marker = fr.path_marker(ca, cb, node_path, curve_path, curves, by_code, nm2)
got = json.loads(marker)
check("PathMarker含卡扣",
      got == ["中控大屏", "Body.46", "C14\\Tyton-157-00181", "Body.47",
              "组合仪表显示屏"], marker)

# ---------- 8. 断路 ----------
nodes2, segs2 = wire(
    [mknode("CA", "connector", 1), mknode("HA", "hanging", 1),
     mknode("CB", "connector", 1), mknode("HB", "hanging", 1)],
    [mkseg("S1", "CA", "HA", 100, [0, 0, 0], [100, 0, 0], "B1"),
     mkseg("S2", "CB", "HB", 100, [500, 0, 0], [600, 0, 0], "B2")])
for n, xyz in zip(nodes2, [[0, 0, 0], [100, 0, 0], [500, 0, 0], [600, 0, 0]]):
    n["xyz"] = [float(v) for v in xyz]
for n in nodes2:
    n["_fx"] = fr.node_face_center(n, {s["seg"]: s for s in segs2})
by2 = {n["code"]: n for n in nodes2}
curves2 = fr.cluster_curves(nodes2, segs2, lambda b: b)
adj2 = fr.build_graph(curves2)
check("跨域不连通", fr.shortest_path(adj2, curves2, "CA", "CB") is None)
ca2 = {"partname": "A", "node": "CA"}
cb2 = {"partname": "B", "node": "CB"}
d = fr.broken_detail(ca2, cb2, adj2, curves2, by2)
check("断路诊断给出gap", "gap=400.0mm" in d, d)
check("断路诊断定位两侧", "HA" in d and "CB" in d, d)

# ---------- 9. 同点次要连接器 ----------
n = mknode("CA", "connector", 1)
n["entity"] = "X#S4"
n["secondary"] = [{"code": "CONX", "kind": "connector", "tag": "X#S9"}]
conns = fr.collect_connectors([n], {"S4": "CONNECTOR\\测试1\\test",
                                    "S9": "CONNECTOR\\测试2\\test"})
check("次要连接器被收录", len(conns) == 2, str(len(conns)))
check("次要连接器同节点", all(c["node"] == "CA" for c in conns))

# ---------- 10. 原始表读取 ----------
import openpyxl
p = "/tmp/test_orig.xlsx"
wb = openpyxl.Workbook()
ws = wb.active
ws.append(["FROM", "TO", "Length"])
ws.append(["中控大屏", "组合仪表显示屏", 1234])
ws.append(["A", "B", None])
wb.save(p)
rows = fr.read_original_table(p)
check("原始表读取", rows == [("中控大屏", "组合仪表显示屏", 1234.0), ("A", "B", None)],
      str(rows))
os.remove(p)

print("\n==== %d PASS, %d FAIL ====" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
