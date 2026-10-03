#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""线束 STEP(AP242) 拓扑分析核心(结构化判定, 不依赖命名).

分支判定完全基于几何结构:
  管状实体 = 有成组等面积平面端盖 + 长径比>=4 + 侧壁面很少
  多实体复合体会被拆成单个实体逐个判定(扁平 STEP 兼容)
  分支 = 管状实体的中线重叠聚类(半壳/分体合并); 中线优先用空间匹配的真线框,
         否则从实体反推; 对不上任何分支的碎线框直接忽略
用法:
  from harness_topology import analyze
  r = analyze(step_path, tol=3.0, out_dir=None, progress=print, force_reverse=False)
命令行:
  python harness_topology.py xxx.stp [容差mm] [--force-reverse]
"""
import re, math, json, os, sys
from collections import defaultdict, Counter
import numpy as np

from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TDocStd import TDocStd_Document
from OCP.TCollection import TCollection_ExtendedString
from OCP.XCAFDoc import XCAFDoc_DocumentTool, XCAFDoc_ShapeTool, XCAFDoc_Location
from OCP.collections import Sequence_TDF_Label
from OCP.TDF import TDF_Label
from OCP.TDataStd import TDataStd_Name
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID
from OCP.TopoDS import TopoDS
TopoDS_Edge_s = TopoDS.Edge
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.TopLoc import TopLoc_Location
from OCP.GeomAbs import GeomAbs_Plane
from OCP.gp import gp_Trsf, gp_Pnt
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex

S = XCAFDoc_ShapeTool
WIRE_MIN_LEN = 20.0   # 可用线框最小总长 mm
ASPECT_MIN = 4.0      # 管状判定最小长径比
# --- 接触式拓扑分析参数(2026-09-27 方案) ---
CONTACT_TOL = 2.0        # 实体接触判定距离 mm
END_TOL = 5.0            # 支撑点距管端 <= 此值视为端部接触 mm
STATION_END_TOL = 5.0    # 搭接投影距主干端头 <= 此值按端点连接处理 mm
STATION_MERGE_TOL = 10.0 # 搭接/固定站位合并: 弧长差与空间距离双条件 mm
TIE_SIZE_MAX = 30.0      # 非管实体最大包围盒尺寸 <= 此值判为扎带/卡扣 mm


def decode_name(s):
    return re.sub(r'\\X2\\([0-9A-Fa-f]+)\\X0\\',
                  lambda m: ''.join(chr(int(m.group(1)[i:i+4], 16))
                                        for i in range(0, len(m.group(1)), 4)), s)


def label_name(lab):
    a = TDataStd_Name()
    return decode_name(a.Get().ToExtString()) if lab.FindAttribute(TDataStd_Name.GetID_s(), a) else "?"


def comp_location(lab):
    la = XCAFDoc_Location()
    if lab.FindAttribute(XCAFDoc_Location.GetID_s(), la):
        return la.Get()
    return TopLoc_Location()


_AUTO_NAMES = ("COMPOUND", "SOLID")


def is_real_name(n):
    return bool(n) and n != "?" and not n.startswith("=>") and n not in _AUTO_NAMES


def count_topo(shp, t):
    e = TopExp_Explorer(shp, t); n = 0
    while e.More(): n += 1; e.Next()
    return n


def edge_length(edge):
    p = GProp_GProps(); BRepGProp.LinearProperties_s(edge, p)
    return p.Mass()


def sample_edge(edge, n=60):
    ad = BRepAdaptor_Curve(TopoDS_Edge_s(edge))
    u0, u1 = ad.FirstParameter(), ad.LastParameter()
    pts = []
    for k in range(n + 1):
        p = ad.Value(u0 + (u1 - u0) * k / n)
        pts.append((p.X(), p.Y(), p.Z()))
    return np.array(pts)


def xform_pts(pts, trsf):
    out = np.empty_like(pts)
    for i, (x, y, z) in enumerate(pts):
        p = gp_Pnt(x, y, z); p.Transform(trsf)
        out[i] = (p.X(), p.Y(), p.Z())
    return out


def solid_volume_bbox(shape, trsf):
    p = GProp_GProps(); BRepGProp.VolumeProperties_s(shape, p)
    box = Bnd_Box(); BRepBndLib.Add_s(shape, box)
    x0, y0, z0 = box.GetXMin(), box.GetYMin(), box.GetZMin()
    x1, y1, z1 = box.GetXMax(), box.GetYMax(), box.GetZMax()
    c = np.array([(x0+x1)/2, (y0+y1)/2, (z0+z1)/2, x1-x0, y1-y0, z1-z0])
    pt = gp_Pnt(c[0], c[1], c[2]); pt.Transform(trsf)
    c[0], c[1], c[2] = pt.X(), pt.Y(), pt.Z()
    return p.Mass(), c


def face_areas(shp):
    out = []
    ex = TopExp_Explorer(shp, TopAbs_FACE)
    while ex.More():
        f = TopoDS.Face(ex.Current())
        ad = BRepAdaptor_Surface(f)
        p = GProp_GProps(); BRepGProp.SurfaceProperties_s(f, p)
        out.append((ad.GetType() == GeomAbs_Plane, p.Mass(), f))
        ex.Next()
    return out


def is_tube_solid(shp):
    """结构化管状判定: 成组等面积平面端盖 + 长径比>=ASPECT_MIN + 侧壁面很少.
    返回 {'dia','length','cap_faces','lat_faces'} 或 None."""
    fa = face_areas(shp)
    vp = GProp_GProps(); BRepGProp.VolumeProperties_s(shp, vp)
    vol = vp.Mass()
    if vol <= 0 or not fa:
        return None
    planar = sorted(a for pl, a, _ in fa if pl)
    lat = [(a, f) for pl, a, f in fa if not pl]
    if len(planar) < 2 or not lat:
        return None
    lat_area = sum(a for a, _ in lat)
    clusters = []  # 平面按面积聚类(2%容差)
    for a in planar:
        for c in clusters:
            if abs(a - c[0]) / c[0] < 0.02:
                c[0] = (c[0]*c[1] + a) / (c[1] + 1); c[1] += 1
                break
        else:
            clusters.append([a, 1])
    cap_groups = [c for c in clusters if c[1] >= 2]
    if not cap_groups:
        return None
    cap_area, n_cap = min(cap_groups, key=lambda c: c[0])
    if n_cap > 8 or len(fa) > 40:
        return None
    if cap_area * n_cap >= 0.3 * lat_area:
        return None
    dia = 2 * math.sqrt(cap_area / math.pi)
    length = vol / cap_area
    if length / dia < ASPECT_MIN:
        return None
    cap_faces = [f for pl, a, f in fa if pl and abs(a - cap_area) / cap_area < 0.02]
    return {"cap_area": cap_area, "n_cap": n_cap, "dia": dia, "length": length,
            "cap_faces": cap_faces, "lat_faces": [f for _, f in lat]}


def circle_center_3d(pts):
    """空间圆拟合求圆心(半圈轮廓也能正确还原圆心); 退化时返回形心."""
    c0 = pts.mean(axis=0)
    X = pts - c0
    try:
        _, _, vt = np.linalg.svd(X, full_matrices=False)
    except Exception:
        return c0
    e1, e2 = vt[0], vt[1]
    p2 = np.column_stack([X @ e1, X @ e2])
    x, y = p2[:, 0], p2[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x**2 + y**2)
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except Exception:
        return c0
    cx, cy = -sol[0]/2, -sol[1]/2
    return c0 + cx*e1 + cy*e2


def profile_circle(pts):
    """轮廓点拟合圆, 返回 (圆心, 半径, 相对残差)."""
    c = circle_center_3d(pts)
    dist = np.linalg.norm(pts - c, axis=1)
    r = float(np.mean(dist))
    res = float(np.std(dist) / max(r, 1e-9))
    return c, r, res


def face_centerline(face, trsf, dia, n_station=48, m_profile=20):
    """侧壁面轴线: 用管径先验判定扫掠方向(轮廓圆拟合质量), 再逐站位圆拟合求圆心."""
    try:
        ad = BRepAdaptor_Surface(face)
        u0, u1 = ad.FirstUParameter(), ad.LastUParameter()
        v0, v1 = ad.FirstVParameter(), ad.LastVParameter()
        if not all(math.isfinite(x) for x in (u0, u1, v0, v1)):
            return None
        umid, vmid = (u0+u1)/2, (v0+v1)/2
        re = max(dia/2, 1e-6)

        def prof_points(along_u):
            pts = []
            for j in range(m_profile):
                s = j/m_profile
                if along_u:  # 扫掠沿U, 轮廓沿V
                    p = ad.Value(umid, v0+(v1-v0)*s)
                else:        # 扫掠沿V, 轮廓沿U
                    p = ad.Value(u0+(u1-u0)*s, vmid)
                g = gp_Pnt(p.X(), p.Y(), p.Z()); g.Transform(trsf)
                pts.append((g.X(), g.Y(), g.Z()))
            return np.array(pts)

        def dir_score(along_u):
            _, r, res = profile_circle(prof_points(along_u))
            if not (math.isfinite(r) and r > 1e-9):
                return 1e9
            return abs(math.log(r/re)) + res

        su, sv = dir_score(True), dir_score(False)
        if min(su, sv) > 5:  # 两边拟合都差, 不是管壁面
            return None
        along_u = su <= sv

        spine = []
        for k in range(n_station+1):
            t = k/n_station
            prof = []
            for j in range(m_profile):
                s = j/m_profile
                if along_u:
                    p = ad.Value(u0+(u1-u0)*t, v0+(v1-v0)*s)
                else:
                    p = ad.Value(u0+(u1-u0)*s, v0+(v1-v0)*t)
                g = gp_Pnt(p.X(), p.Y(), p.Z()); g.Transform(trsf)
                prof.append((g.X(), g.Y(), g.Z()))
            prof = np.array(prof)
            c, r, res = profile_circle(prof)
            if not np.all(np.isfinite(c)) or not (0.3*re < r < 3*re) or res > 0.5:
                c = prof.mean(axis=0)  # 退化回形心
            spine.append(c)
        spine = np.array(spine)
        if not np.all(np.isfinite(spine)):
            return None
        return spine
    except Exception:
        return None


def _end_tangent(pts, end, window=5.0):
    """端点处切向(向外为正), 取端点附近 window(mm)内的点拟合方向."""
    if len(pts) < 2:
        return None
    p = pts if end == 1 else pts[::-1]
    tip = p[-1]
    d = np.linalg.norm(p - tip, axis=1)
    sel = p[d <= window]
    if len(sel) < 2:
        sel = p[-2:]
    v = sel[-1] - sel[0]
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else None


def chain_all(curves, tol=1.0, max_angle_deg=None):
    """把零散曲线段按端点拼接成链, 返回全部链 [{"pts","len","n_seg","members"}] (按长度降序).
    扁平 STEP 里一个产品可能含多条互不相连的线框, 需保留每条链分别判定.
    max_angle_deg: 拼接时两段在接头处的切向夹角上限(度); 用于区分"同一扫掠体的分段"
    (相切, 应合并)与"分叉节点的两根管子"(有夹角, 不合并). None 表示不检查角度."""
    segs = [{"pts": c["pts"], "len": c["len"], "src": c} for c in curves]
    if not segs:
        return []
    max_cos = math.cos(math.radians(max_angle_deg)) if max_angle_deg else None
    used = [False]*len(segs)
    order = sorted(range(len(segs)), key=lambda i: -segs[i]["len"])
    out = []
    for s0 in order:
        if used[s0]:
            continue
        used[s0] = True
        pts = segs[s0]["pts"]
        members = [segs[s0]["src"]]
        nseg = 1
        for end in (0, 1):
            while True:
                tip = pts[0] if end == 0 else pts[-1]
                best, bestd, flip = -1, tol, False
                for j, s in enumerate(segs):
                    if used[j]:
                        continue
                    q = s["pts"]
                    d0 = np.linalg.norm(q[0]-tip)
                    d1 = np.linalg.norm(q[-1]-tip)
                    bd = d0 if d0 <= d1 else d1
                    if bd >= bestd:
                        continue
                    # flip 使 q 的"连接端"朝向 tip, 另一端向外延伸
                    # end=0(前端): 连接端应为 q[-1] -> d0近则翻转
                    # end=1(后端): 连接端应为 q[0]  -> d1近则翻转
                    fj = (d0 <= d1) if end == 0 else (d1 < d0)
                    if max_cos is not None:
                        t_tip = _end_tangent(pts, end)
                        # q 的连接端: end=1 时不翻转接 q[0]、翻转接 q[-1];
                        # end=0 时翻转接 q[0]、不翻转接 q[-1]
                        ce = (0 if not fj else 1) if end == 1 else (0 if fj else 1)
                        t_qe = _end_tangent(q, ce)  # 连接端处"向外"(指向接头)切向
                        # 相切延续要求: 链的向外切向 t_tip 与 q 从接头出发的行进方向
                        # (-t_qe) 夹角小, 即 t_tip·t_qe <= -cos(上限); 否则是分叉, 不拼接
                        if t_tip is not None and t_qe is not None:
                            if float(t_tip @ t_qe) > -max_cos:
                                continue
                    flip, best, bestd = fj, j, bd
                if best < 0:
                    break
                used[best] = True
                nseg += 1
                members.append(segs[best]["src"])
                q = segs[best]["pts"][::-1] if flip else segs[best]["pts"]
                pts = np.vstack([q[:-1], pts]) if end == 0 else np.vstack([pts, q[1:]])
        L = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        out.append({"pts": pts, "len": L, "n_seg": nseg, "members": members})
    out.sort(key=lambda d: -d["len"])
    return out


def chain_curves(curves, tol=1.0):
    chains = chain_all(curves, tol)
    if not chains:
        return np.zeros((0, 3)), 0.0
    return chains[0]["pts"], chains[0]["len"]


def polyline_dist(pts, poly, n=60):
    """pts 各点到折线 poly 的最短距离."""
    if len(poly) < 2 or len(pts) == 0:
        return np.full(len(pts), 1e9)
    a = poly[:-1]; b = poly[1:]
    ab = b - a
    denom = np.einsum('ij,ij->i', ab, ab) + 1e-12
    out = np.empty(len(pts))
    for i, p in enumerate(pts):
        t = np.clip(np.einsum('ij,ij->i', ab, p - a) / denom, 0, 1)
        proj = a + t[:, None]*ab
        out[i] = np.linalg.norm(proj - p, axis=1).min()
    return out


def spines_overlap(a, b, tol=2.0, frac=0.6):
    """两条中线是否大面积重叠(半壳/重复分面)."""
    if len(a) < 2 or len(b) < 2:
        return False
    short, long = (a, b) if len(a) <= len(b) else (b, a)
    d = polyline_dist(short[::max(1, len(short)//60)], long)
    return np.mean(d < tol) >= frac


def reverse_spine(tube, trsf):
    """从管状实体反推中心线: 各侧壁参数中线 -> 拼接 -> 端点对齐端盖中心."""
    segs = []
    for f in tube["lat_faces"]:
        sp = face_centerline(f, trsf, tube["dia"])
        if sp is not None and len(sp) > 1 and np.all(np.isfinite(sp)):
            L = float(np.sum(np.linalg.norm(np.diff(sp, axis=0), axis=1)))
            if L > 1e-6:
                segs.append({"pts": sp, "len": L})
    if not segs:
        return None
    # 去重: 半壳等重复分面的中线相互重叠, 只保留最长的
    kept = []
    for s in sorted(segs, key=lambda s: -s["len"]):
        if any(spines_overlap(s["pts"], k["pts"]) for k in kept):
            continue
        kept.append(s)
    pts, _ = chain_curves(kept, tol=2.0)
    caps = []
    for f in tube["cap_faces"]:
        p = GProp_GProps(); BRepGProp.SurfaceProperties_s(f, p)
        c = p.CentreOfMass(); g = gp_Pnt(c.X(), c.Y(), c.Z()); g.Transform(trsf)
        caps.append(np.array([g.X(), g.Y(), g.Z()]))
    if len(pts) >= 2 and caps:
        for end in (0, -1):
            d = [np.linalg.norm(c - pts[end]) for c in caps]
            j = int(np.argmin(d))
            if d[j] < 15.0:
                pts = np.vstack([caps[j], pts[1:]]) if end == 0 else np.vstack([pts[:-1], caps[j]])
    L = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    return pts, L


# ======================================================================
# 接触式拓扑分析(2026-09-27 方案):
#   管-管接触 -> 端-端连接 / T 型搭接(投影站位)
#   管-连接器/扎带接触 -> 端部命名节点 / 中部固定站位
#   站位聚类(弧长+空间双条件 10mm) -> 切分线段 -> 节点 -> 全局编码 -> 走线
# ======================================================================

def _d3(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def _poly_project(p, pts):
    """点 p 投影到折线 pts. 返回 (弧长站位 s, 投影点 xyz, 距离)."""
    best = None
    s = 0.0
    for i in range(len(pts) - 1):
        ax, ay, az = pts[i]
        cx, cy, cz = pts[i+1]
        vx, vy, vz = cx-ax, cy-ay, cz-az
        L2 = vx*vx + vy*vy + vz*vz
        t = 0.0 if L2 < 1e-12 else ((p[0]-ax)*vx + (p[1]-ay)*vy + (p[2]-az)*vz) / L2
        t = max(0.0, min(1.0, t))
        q = (ax+vx*t, ay+vy*t, az+vz*t)
        d = _d3(p, q)
        seglen = math.sqrt(L2)
        if best is None or d < best[2]:
            best = (s + seglen*t, q, d)
        s += seglen
    return best


def _spine_point_at(pts, s):
    """折线上弧长 s 处的点(线性插值)."""
    acc = 0.0
    for i in range(len(pts)-1):
        a, b = pts[i], pts[i+1]
        L = _d3(a, b)
        if acc + L >= s or i == len(pts)-2:
            t = 0.0 if L < 1e-12 else max(0.0, min(1.0, (s-acc)/L))
            return (a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t, a[2]+(b[2]-a[2])*t)
        acc += L
    return tuple(pts[-1])


def _cluster_stations(stations, tol):
    """站位聚类: |弧长差|<=tol 且 空间距离<=tol 才合并(双条件).
    stations: [{'s','xyz','kind',...}]. 返回 [{'s','xyz','members'}]."""
    groups = []
    for st in sorted(stations, key=lambda r: r["s"]):
        hit = None
        for g in groups:
            if abs(g["s"] - st["s"]) <= tol and _d3(g["xyz"], st["xyz"]) <= tol:
                hit = g
                break
        if hit is None:
            groups.append({"s": st["s"], "xyz": tuple(st["xyz"]), "members": [st]})
        else:
            m = hit["members"] + [st]
            hit["members"] = m
            hit["s"] = sum(r["s"] for r in m) / len(m)
            hit["xyz"] = tuple(sum(r["xyz"][k] for r in m)/len(m) for k in range(3))
    return groups


def _split_polyline(pts, stations):
    """在站位(弧长 s 列表)处切分折线. 返回 [{'pts','s0','s1'}]."""
    total = sum(_d3(pts[i], pts[i+1]) for i in range(len(pts)-1))
    ss = sorted(set(round(s, 3) for s in stations))
    cuts = [0.0] + [s for s in ss if 1e-6 < s < total - 1e-6] + [total]
    cum = [0.0]
    for i in range(len(pts)-1):
        cum.append(cum[-1] + _d3(pts[i], pts[i+1]))
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a < 1e-6:
            continue
        seg = [_spine_point_at(pts, a)]
        for i, c in enumerate(cum):
            if a + 1e-9 < c < b - 1e-9:
                seg.append(tuple(pts[i]))
        seg.append(_spine_point_at(pts, b))
        out.append({"pts": seg, "s0": a, "s1": b})
    return out


def _world_shape(shape, trsf):
    """返回施加变换后的形状(不污染被共享的原型)."""
    from OCP.TopLoc import TopLoc_Location
    return shape.Located(TopLoc_Location(trsf))


def _shape_contact(sa, sb):
    """两实体最小距离与支撑点. 返回 (dist, pa_xyz, pb_xyz) 或 None."""
    dss = BRepExtrema_DistShapeShape(sa, sb)
    if not dss.Perform():
        return None
    if not dss.IsDone() or dss.NbSolution() < 1:
        return None
    p1, p2 = dss.PointOnShape1(1), dss.PointOnShape2(1)
    return (dss.Value(), (p1.X(), p1.Y(), p1.Z()), (p2.X(), p2.Y(), p2.Z()))


def _end_touch_other(end_pt, other_ws):
    """管端面中心到另一分支各实体的最小距离(点到实体).
    对互穿/嵌入式接触鲁棒: 嵌入的端面中心落在对方实体内(距离0).
    返回最小距离, 失败返回 None."""
    best = None
    try:
        v = BRepBuilderAPI_MakeVertex(gp_Pnt(float(end_pt[0]), float(end_pt[1]),
                                             float(end_pt[2]))).Vertex()
    except Exception:
        return None
    for b in other_ws:
        try:
            dss = BRepExtrema_DistShapeShape(v, b["w"])
            if dss.Perform() and dss.IsDone() and dss.NbSolution() >= 1:
                d = dss.Value()
                if best is None or d < best:
                    best = d
        except Exception:
            continue
    return best


def _spine_joint(spA, spB, tol=2.0):
    """两条 spine 的端点若在 tol 内相接(拼接规则), 返回接头中点, 否则 None."""
    for ea in (0, 1):
        pa = tuple(spA[0] if ea == 0 else spA[-1])
        for eb in (0, 1):
            pb = tuple(spB[0] if eb == 0 else spB[-1])
            if _d3(pa, pb) <= tol:
                return ((pa[0]+pb[0])/2, (pa[1]+pb[1])/2, (pa[2]+pb[2])/2)
    return None


def _item_radius(it):
    v = it["solid"]["vol"]; L = it["spine_len"]
    return math.sqrt(v/(math.pi*L)) if v and L and L > 0 else 5.0


def _merge_entity_tag(tags):
    """共位实体合并 tag: ['P#S7','P#S8'] -> 'P#S7+S8'."""
    pres, nums, plain = set(), [], []
    for t in tags:
        if "#S" in t:
            pre, num = t.rsplit("#S", 1)
            pres.add(pre); nums.append(num)
        else:
            plain.append(t)
    if pres and not plain and len(pres) == 1:
        nums = sorted(set(nums), key=lambda x: int(x) if x.isdigit() else x)
        return pres.pop() + "#S" + "+S".join(nums)
    return "+".join(sorted(tags))


def _alias_colocated(conn_candidates, log):
    """共位非管实体合并: 包围盒中心 <=10mm 视为同一逻辑件(如同一管端的两片式).
    返回 {原tag: 合并后tag}."""
    tags = [s["tag"] for s in conn_candidates]
    centers = {s["tag"]: tuple(s["bbox_c"][:3]) for s in conn_candidates}
    alias, used = {}, set()
    for t in sorted(tags):
        if t in used:
            continue
        grp = [t]; used.add(t)
        for u in sorted(tags):
            if u in used:
                continue
            if _d3(centers[t], centers[u]) <= 10.0:
                grp.append(u); used.add(u)
        new_tag = _merge_entity_tag(grp) if len(grp) > 1 else grp[0]
        for g in grp:
            alias[g] = new_tag
        if len(grp) > 1:
            log(f"共位实体合并: {' + '.join(sorted(grp))} -> {new_tag}")
    return alias


def _classify_entities(branch_data, relations, log):
    """非管实体拓扑分类(用户规则):
       连接器 = 只接触一个线段端头;
       固定卡扣 = 线束从中穿过(中部站位)或接触多个端头.
       端-端接头展开: 实体若是接头处唯一的设备, 则把接头另一侧端头也计入
       (线束路径从中穿过 -> 路径中间 -> 固定卡扣);
       若接头处还有别的设备, 只计直接接触的端头."""
    parent = {}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for (i, ei, j, ej) in relations.get("end_ends", []):
        for kk in ((i, ei), (j, ej)):
            parent.setdefault(kk, kk)
        ra, rb = find((i, ei)), find((j, ej))
        if ra != rb:
            parent[rb] = ra
    end_tags = defaultdict(set)  # (branch,end) -> 直接接触的实体tag集合
    for t in relations.get("terminals", []):
        end_tags[(t["branch"], t["end"])].add(t["tag"])
    term_ends = defaultdict(set)
    for t in relations.get("terminals", []):
        be = (t["branch"], t["end"])
        term_ends[t["tag"]].add(be)
        if be in parent:
            r = find(be)
            others = [k for k in parent if find(k) == r and k != be]
            if others and not any(end_tags[o] - {t["tag"]} for o in others):
                term_ends[t["tag"]].update(others)
    stat_tags = set(t["tag"] for t in relations.get("tie_stations", []))
    kind = {}
    for tag in set(term_ends) | stat_tags:
        n_end = len(term_ends.get(tag, ()))
        if tag in stat_tags or n_end != 1:
            kind[tag] = "clamp"
            if n_end > 1:
                log(f"分类: {tag} 接触{n_end}个线段端头 -> 固定卡扣")
        else:
            kind[tag] = "connector"
    return kind


def compute_relations(branch_data, conn_candidates, log):
    """实体接触分析(需 OCP).
    返回 {"taps":[{main,s,xyz,tap,tap_end,dev}],
            "end_ends":[(bi,ei,bj,ej)],
            "terminals":[{branch,end,tag,kind,dist}],
            "tie_stations":[{branch,s,xyz,tag,kind,dev}],
            "side_sides":[{a,b,dist}]}."""
    rel = {"taps": [], "end_ends": [], "terminals": [], "tie_stations": [],
           "side_sides": []}

    def _world_box(w):
        """世界坐标系下的包围盒 (用于 BRepExtrema 前的精确粗筛)."""
        box = Bnd_Box()
        BRepBndLib.Add_s(w, box)
        return box

    bw = []
    for b in branch_data:
        ws = []
        for sd in b.get("solids", []):
            try:
                w = _world_shape(sd["shape"], sd["trsf"])
            except Exception:
                continue
            c = sd["bbox_c"]
            ws.append({"w": w, "c": tuple(c[:3]),
                       "r": float(np.linalg.norm(c[3:]))/2.0,
                       "box": _world_box(w)})
        bw.append(ws)
    cw = []
    for s in conn_candidates:
        try:
            w = _world_shape(s["shape"], s["trsf"])
        except Exception:
            continue
        c = s["bbox_c"]
        cw.append({"w": w, "tag": s["tag"], "proto": s["proto"],
                   "c": tuple(c[:3]), "r": float(np.linalg.norm(c[3:]))/2.0,
                   "box": _world_box(w),
                   "kind": "entity"})  # 实体类别由拓扑分类确定, 此处仅占位

    def broad(a, b):
        return _d3(a["c"], b["c"]) <= a["r"] + b["r"] + CONTACT_TOL

    def min_contact(la, lb):
        best = None
        for a in la:
            for b in lb:
                if not broad(a, b):
                    continue
                # 包围盒精确距离粗筛: 盒距 > CONTACT_TOL 的对精确距离必 > 阈值,
                # 本来也会被调用方跳过, 直接跳过 BRepExtrema (结果完全一致).
                if a["box"].Distance(b["box"]) > CONTACT_TOL:
                    continue
                r = _shape_contact(a["w"], b["w"])
                if r and (best is None or r[0] < best[0]):
                    best = r
        return best

    def spine_len(pts):
        return sum(_d3(pts[i], pts[i+1]) for i in range(len(pts)-1))

    n = len(branch_data)
    branch_radius = []
    for b in branch_data:
        v = b.get("vol"); L = b.get("length")
        branch_radius.append(math.sqrt(v/(math.pi*L)) if v and L and L > 0 else 5.0)
    jtap_seen = set()
    # 管-管
    for i in range(n):
        if not bw[i]:
            continue
        pi = branch_data[i]
        ptsi = [tuple(p) for p in pi["pts"]]
        for j in range(i+1, n):
            if not bw[j]:
                continue
            r = min_contact(bw[i], bw[j])
            if not r or r[0] > CONTACT_TOL:
                continue
            dist = r[0]
            pj = branch_data[j]
            ptsj = [tuple(p) for p in pj["pts"]]
            # 端面中心接触判定: 对互穿/嵌入式搭接鲁棒(支撑点法在嵌入时会失效)
            hi = [(e, d) for e in (0, 1)
                  for d in [_end_touch_other(ptsi[0] if e == 0 else ptsi[-1], bw[j])]
                  if d is not None and d <= CONTACT_TOL]
            hj = [(e, d) for e in (0, 1)
                  for d in [_end_touch_other(ptsj[0] if e == 0 else ptsj[-1], bw[i])]
                  if d is not None and d <= CONTACT_TOL]
            if hi and hj:
                ei = min(hi, key=lambda t: t[1])[0]
                ej = min(hj, key=lambda t: t[1])[0]
                rel["end_ends"].append((i, ei, j, ej))
                log(f"接触: {pi['key']}端 <-> {pj['key']}端 ({dist:.1f}mm)")
                # 接头若落在第三条管体中部 -> 双分支搭接点(degree=4 四路节点):
                # 端-端连接保留, 同时在主干上记一处搭接站位
                pe_i = ptsi[0] if ei == 0 else ptsi[-1]
                pe_j = ptsj[0] if ej == 0 else ptsj[-1]
                J = ((pe_i[0]+pe_j[0])/2, (pe_i[1]+pe_j[1])/2, (pe_i[2]+pe_j[2])/2)
                for k in range(n):
                    if k == i or k == j or not bw[k]:
                        continue
                    pk = branch_data[k]
                    ptsk = [tuple(p) for p in pk["pts"]]
                    s, q, dev = _poly_project(J, ptsk)
                    Lk = spine_len(ptsk)
                    if (dev <= branch_radius[i] + branch_radius[k] + CONTACT_TOL
                            and STATION_END_TOL < s < Lk - STATION_END_TOL
                            and (k, i, j) not in jtap_seen):
                        jtap_seen.add((k, i, j))
                        rel["taps"].append({"main": k, "s": s, "xyz": q, "tap": i,
                                            "tap_end": ei, "dev": dev})
                        log(f"接头搭接(四路节点): {pi['key']}+{pj['key']}接头 -> "
                            f"{pk['key']} 站位{s:.1f}mm 偏差{dev:.1f}mm")
                        break
            elif hi:
                ei = min(hi, key=lambda t: t[1])[0]
                pe = ptsi[0] if ei == 0 else ptsi[-1]
                s, q, dev = _poly_project(pe, ptsj)
                Lj = spine_len(ptsj)
                if s <= STATION_END_TOL or s >= Lj - STATION_END_TOL:
                    rel["end_ends"].append((i, ei, j, 0 if s < Lj/2 else 1))
                    log(f"接触(近端按端点连接): {pi['key']} <-> {pj['key']}端 ({dist:.1f}mm)")
                else:
                    rel["taps"].append({"main": j, "s": s, "xyz": q, "tap": i,
                                        "tap_end": ei, "dev": dev})
                    log(f"搭接: {pi['key']} -> {pj['key']} 站位{s:.1f}mm 偏差{dev:.1f}mm")
            elif hj:
                ej = min(hj, key=lambda t: t[1])[0]
                pe = ptsj[0] if ej == 0 else ptsj[-1]
                s, q, dev = _poly_project(pe, ptsi)
                Li = spine_len(ptsi)
                if s <= STATION_END_TOL or s >= Li - STATION_END_TOL:
                    rel["end_ends"].append((j, ej, i, 0 if s < Li/2 else 1))
                    log(f"接触(近端按端点连接): {pj['key']} <-> {pi['key']}端 ({dist:.1f}mm)")
                else:
                    rel["taps"].append({"main": i, "s": s, "xyz": q, "tap": j,
                                        "tap_end": ej, "dev": dev})
                    log(f"搭接: {pj['key']} -> {pi['key']} 站位{s:.1f}mm 偏差{dev:.1f}mm")
            else:
                rel["side_sides"].append({"a": pi["key"], "b": pj["key"], "dist": dist})
                log(f"侧-侧接触(忽略, 手工处理): {pi['key']} <-> {pj['key']} ({dist:.1f}mm)")
    # 管-非管(连接器/扎带): 两遍扫描, 先端部后中部.
    # 同一刚体在端部已命名后, 其在附近管段上的附带接触不再另起站位
    # (避免同一连接器/扎带被编出两个 CON/TIE 码, 破坏编码对照表 1:1).
    term_xyz = []  # (tag, 端部xyz, 包围盒对角线)
    for i in range(n):
        if not bw[i]:
            continue
        pi = branch_data[i]
        ptsi = [tuple(p) for p in pi["pts"]]
        Li = spine_len(ptsi)
        for c in cw:
            r = min_contact(bw[i], [c])
            if not r or r[0] > CONTACT_TOL:
                continue
            dist, pa, _pc = r
            s, q, dev = _poly_project(pa, ptsi)
            if s <= END_TOL or s >= Li - END_TOL:
                end = 0 if s < Li/2 else 1
                exyz = ptsi[0] if end == 0 else ptsi[-1]
                rel["terminals"].append({"branch": i, "end": end, "tag": c["tag"],
                                        "kind": c["kind"], "dist": dist})
                term_xyz.append((c["tag"], exyz, 2.0*c["r"]))
                log(f"端部接触: {pi['key']}[{end}] <-> 实体 {c['tag']} ({dist:.1f}mm)")
    for i in range(n):
        if not bw[i]:
            continue
        pi = branch_data[i]
        ptsi = [tuple(p) for p in pi["pts"]]
        Li = spine_len(ptsi)
        for c in cw:
            r = min_contact(bw[i], [c])
            if not r or r[0] > CONTACT_TOL:
                continue
            dist, pa, _pc = r
            s, q, dev = _poly_project(pa, ptsi)
            if s <= END_TOL or s >= Li - END_TOL:
                continue  # 第一遍已处理
            dup = next((t for t in term_xyz if t[0] == c["tag"]
                        and _d3(q, t[1]) <= t[2]), None)
            if dup is not None:
                log(f"同体接触(跳过): {pi['key']} 站位{s:.1f}mm <-> 实体 {c['tag']} "
                    f"(距已命名端部{_d3(q, dup[1]):.1f}mm, 同一刚体)")
                continue
            rel["tie_stations"].append({"branch": i, "s": s, "xyz": q,
                                       "tag": c["tag"], "kind": c["kind"], "dev": dev})
            log(f"中部固定: {pi['key']} 站位{s:.1f}mm <-> 实体 {c['tag']} ({dist:.1f}mm)")
    return rel


def build_topology_from_relations(branch_data, relations, tol, conn_candidates, log):
    """由接触关系构建拓扑(纯逻辑, 无 OCP):
       站位聚类 -> 切分线段 -> 节点(含强制归属与传递合并) -> 全局编码 -> 走线.
    返回 (segments, nodes, entities, runs, branch_info)."""
    n = len(branch_data)
    # 0) 共位实体合并 + 非管实体拓扑分类(用户规则: 连接器只连一个线段端头)
    alias = _alias_colocated(conn_candidates, log)
    for t in relations.get("terminals", []):
        t["tag"] = alias.get(t["tag"], t["tag"])
    for t in relations.get("tie_stations", []):
        t["tag"] = alias.get(t["tag"], t["tag"])
    kind_of = _classify_entities(branch_data, relations, log)
    for t in relations.get("terminals", []):
        if t["tag"] in kind_of:
            t["kind"] = kind_of[t["tag"]]
    for t in relations.get("tie_stations", []):
        if t["tag"] in kind_of:
            t["kind"] = kind_of[t["tag"]]
    relations["_tag_alias"] = {k: v for k, v in alias.items() if k != v}
    relations["_entity_kind"] = kind_of
    # 1) 每分支站位收集
    stations = [[] for _ in range(n)]
    for t in relations.get("taps", []):
        stations[t["main"]].append({"s": t["s"], "xyz": tuple(t["xyz"]), "kind": "tap",
                                    "tap": t["tap"], "tap_end": t["tap_end"]})
    for t in relations.get("tie_stations", []):
        stations[t["branch"]].append({"s": t["s"], "xyz": tuple(t["xyz"]), "kind": "tie",
                                      "tag": t["tag"], "ekind": t["kind"]})
    # 2) 站位聚类(弧长+空间双条件) -> 标准站位(投影点重算, 保证落在中心线上)
    std_stations = []
    for bi in range(n):
        pts = [tuple(p) for p in branch_data[bi]["pts"]]
        total = sum(_d3(pts[i], pts[i+1]) for i in range(len(pts)-1))
        std = []
        for g in _cluster_stations(stations[bi], STATION_MERGE_TOL):
            s = max(0.0, min(total, g["s"]))
            std.append({"s": s, "xyz": _spine_point_at(pts, s), "members": g["members"],
                        "interior": 1e-6 < s < total - 1e-6})
        std_stations.append(sorted(std, key=lambda r: r["s"]))
    log(f"搭接/固定站位: {sum(len(s) for s in std_stations)}")
    # 3) 切分 -> 线段
    segments = []
    branch_seg = [[] for _ in range(n)]
    interior_of = []
    for bi in range(n):
        pts = [tuple(p) for p in branch_data[bi]["pts"]]
        interior = [st for st in std_stations[bi] if st["interior"]]
        interior_of.append(interior)
        for part in _split_polyline(pts, [st["s"] for st in interior]):
            sid = len(segments)
            segments.append({"id": sid, "branch": bi, "pts": part["pts"],
                             "length": part["s1"]-part["s0"],
                             "p0": part["pts"][0], "p1": part["pts"][-1]})
            branch_seg[bi].append(sid)
    # 4) 节点
    nodes = []  # {xyz, fixed, members:[(seg,end)], entities:[(kind,tag)]}
    force = defaultdict(list)  # (seg,end) -> [node idx] 接触关系强制归属
    assigned = set()

    def new_node(xyz, fixed):
        nodes.append({"xyz": xyz, "fixed": fixed, "members": [], "entities": []})
        return len(nodes)-1

    for bi, std in enumerate(std_stations):
        int_st = interior_of[bi]
        for k, st in enumerate(int_st):
            ni = new_node(st["xyz"], True)
            ka, kb = (branch_seg[bi][k], 1), (branch_seg[bi][k+1], 0)
            nodes[ni]["members"].extend([ka, kb])
            assigned.update([ka, kb])
            for m in st["members"]:
                if m["kind"] == "tie":
                    nodes[ni]["entities"].append((m["ekind"], m["tag"]))
                else:
                    ti, te = m["tap"], m["tap_end"]
                    key = (branch_seg[ti][0] if te == 0 else branch_seg[ti][-1], te)
                    force[key].append(ni)
        for st in std:
            if st["interior"]:
                continue
            ni = new_node(st["xyz"], True)  # 端头站位: 后续与端点节点合并
            for m in st["members"]:
                if m["kind"] == "tie":
                    nodes[ni]["entities"].append((m["ekind"], m["tag"]))
    # 端-端接触组(传递闭包)
    parent = {}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for (i, ei, j, ej) in relations.get("end_ends", []):
        for kk in ((i, ei), (j, ej)):
            parent.setdefault(kk, kk)
        ra, rb = find((i, ei)), find((j, ej))
        if ra != rb:
            parent[rb] = ra
    eegroup = defaultdict(list)
    for kk in parent:
        eegroup[find(kk)].append(kk)
    for grp in eegroup.values():
        ni = new_node(None, False)
        for (bi, e) in grp:
            key = (branch_seg[bi][0] if e == 0 else branch_seg[bi][-1], e)
            force[key].append(ni)
            nodes[ni]["members"].append(key)
            assigned.add(key)
    for nd in nodes:
        if nd["xyz"] is None:
            m = nd["members"]
            nd["xyz"] = tuple(sum(segments[s]["p0" if e == 0 else "p1"][k]
                                  for s, e in m)/len(m) for k in range(3))
    # 其余线段端点: 强制归属优先, 否则空间聚类
    for sg in segments:
        for e in (0, 1):
            key = (sg["id"], e)
            if key in assigned:
                continue
            p = sg["p0"] if e == 0 else sg["p1"]
            if key in force:
                for ni in force[key]:
                    nodes[ni]["members"].append(key)
                assigned.add(key)
                continue
            hit = None
            for ni, nd in enumerate(nodes):
                if _d3(nd["xyz"], p) <= tol:
                    hit = ni
                    break
            if hit is None:
                hit = new_node(tuple(p), False)
            nodes[hit]["members"].append(key)
            assigned.add(key)
    # 节点合并: D1 同一(seg,end)被强制到多节点 -> 同一物理点; D2 空间容差传递闭包
    npar = list(range(len(nodes)))
    def nfind(a):
        while npar[a] != a:
            npar[a] = npar[npar[a]]
            a = npar[a]
        return a
    def nunion(a, b):
        ra, rb = nfind(a), nfind(b)
        if ra != rb:
            npar[rb] = ra
    for key, nilist in force.items():
        for ni in nilist[1:]:
            nunion(nilist[0], ni)
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            if _d3(nodes[i]["xyz"], nodes[j]["xyz"]) <= tol:
                nunion(i, j)
    # D3: 同一连接器/扎带实体只对应一个节点(刚体去重, 保编码对照表 1:1)
    ediag = {}
    for s in conn_candidates:
        try:
            ediag[s["tag"]] = float(np.linalg.norm(s["bbox_c"][3:]))
        except Exception:
            pass
    ent_nodes = defaultdict(list)
    for ni, nd in enumerate(nodes):
        for k, tag in nd["entities"]:
            ent_nodes[(k, tag)].append(ni)
    for (k, tag), nilist in ent_nodes.items():
        if len(nilist) < 2:
            continue
        diag = ediag.get(tag, 15.0)
        for a in range(len(nilist)):
            for b in range(a + 1, len(nilist)):
                ia, ib = nilist[a], nilist[b]
                if nfind(ia) == nfind(ib):
                    continue
                if _d3(nodes[ia]["xyz"], nodes[ib]["xyz"]) <= diag:
                    nunion(ia, ib)
                    log(f"同体节点合并: {tag} ({_d3(nodes[ia]['xyz'], nodes[ib]['xyz']):.1f}mm)")
    merged, root_map = [], {}
    for i, nd in enumerate(nodes):
        r = nfind(i)
        if r not in root_map:
            root_map[r] = len(merged)
            merged.append({"xyz": None, "fixed": False, "members": [], "entities": []})
        m = merged[root_map[r]]
        m["members"].extend(nd["members"])
        m["entities"].extend(nd["entities"])
        if nd["fixed"]:
            m["xyz"] = nd["xyz"]
            m["fixed"] = True
    nodes = merged
    # 无成员节点(如端头站位未合并): 实体信息迁移到最近节点
    keep = []
    for ni, nd in enumerate(nodes):
        if nd["members"]:
            keep.append(nd)
            continue
        best, bestd = None, 1e18
        for md in nodes:
            if not md["members"]:
                continue
            d = _d3(nd["xyz"], md["xyz"])
            if d < bestd:
                bestd, best = d, md
        if best is not None and bestd <= 15.0:
            best["entities"].extend(nd["entities"])
            log(f"孤立站位节点并入 {bestd:.1f}mm 外节点")
        else:
            log(f"警告: 孤立节点无成员, 坐标={nd['xyz']}")
            keep.append(nd)
    nodes = keep
    seg_end_node = {}
    for ni, nd in enumerate(nodes):
        for key in nd["members"]:
            seg_end_node[key] = ni
    # 5) 端部实体关联 -> 命名
    for t in relations.get("terminals", []):
        bi, e = t["branch"], t["end"]
        key = (branch_seg[bi][0] if e == 0 else branch_seg[bi][-1], e)
        ni = seg_end_node.get(key)
        if ni is not None:
            nodes[ni]["entities"].append((t["kind"], t["tag"]))
    counters = {"CON": 0, "TIE": 0, "CLP": 0, "BN": 0, "N": 0}
    for ni, nd in enumerate(nodes):
        seen, um = set(), []
        for key in nd["members"]:
            if key not in seen:
                seen.add(key)
                um.append(key)
        nd["members"] = um
        deg = len(um)
        kinds = {}
        for k, tag in nd["entities"]:
            kinds.setdefault(k, tag)
        # 同点多设备: 按度数定主从. 度数>=2 的节点由穿过型设备(卡扣)主导,
        # 连接器只连一个线段端头, 不能主导多线段节点(用户规则).
        order = (["clamp", "connector", "tie"] if deg >= 2
                 else ["connector", "clamp", "tie"])
        primary = next((k for k in order if k in kinds), None)
        secondary = [(k, tag) for k, tag in kinds.items() if k != primary]
        if primary == "connector":
            counters["CON"] += 1
            code, ntype = f"CON{counters['CON']:02d}", "connector"
            name_3d, entity = kinds["connector"], kinds["connector"]
        elif primary == "clamp":
            counters["CLP"] += 1
            code, ntype = f"CLP{counters['CLP']:02d}", "clamp"
            name_3d, entity = kinds["clamp"], kinds["clamp"]
        elif primary == "tie":
            counters["TIE"] += 1
            code, ntype = f"TIE{counters['TIE']:02d}", "tie"
            name_3d, entity = kinds["tie"], kinds["tie"]
        elif deg >= 3:
            counters["BN"] += 1
            code, ntype = f"BN{counters['BN']:02d}", "fork"
            name_3d, entity = "", None
        else:
            counters["N"] += 1
            code, ntype = f"N{counters['N']:02d}", ("hanging" if deg == 1 else "joint")
            name_3d, entity = "", None
        nd.update({"id": ni, "code": code, "type": ntype, "name_3d": name_3d,
                   "entity": entity, "degree": deg, "secondary": secondary})
        # 次要设备单独编码(同点, 非拓扑节点)
        sec_codes = []
        for sk, stag in secondary:
            if sk == "connector":
                counters["CON"] += 1
                scode = f"CON{counters['CON']:02d}"
            elif sk == "clamp":
                counters["CLP"] += 1
                scode = f"CLP{counters['CLP']:02d}"
            else:
                counters["TIE"] += 1
                scode = f"TIE{counters['TIE']:02d}"
            sec_codes.append((scode, sk, stag))
            log(f"同点多设备: {stag} -> {scode}(次要, 与{code}同位置)")
        nd["secondary_codes"] = sec_codes
        if not nd["fixed"]:
            nd["xyz"] = tuple(sum(segments[s]["p0" if e == 0 else "p1"][k]
                                  for s, e in um)/len(um) for k in range(3)) if um else nd["xyz"]
        nd["xyz"] = tuple(float(x) for x in nd["xyz"])
        nd["nearest"] = None
        if ntype == "hanging":
            best, bestd = None, 1e18
            for s in conn_candidates:
                d = _d3(s["bbox_c"][:3], nd["xyz"])
                sd = max(0.0, d - float(np.linalg.norm(s["bbox_c"][3:]))/2)
                if sd < bestd:
                    bestd, best = sd, s
            if best is not None and bestd <= 30.0:
                nd["nearest"] = {"tag": best["tag"], "proto": best["proto"],
                                 "dist": round(bestd, 1)}
    # 6) 线段端点 -> 节点编码; 走线(度数!=2 处截断); 全局 SEG 编码
    for sg in segments:
        sg["na"] = nodes[seg_end_node[(sg["id"], 0)]]["code"]
        sg["nb"] = nodes[seg_end_node[(sg["id"], 1)]]["code"]
    code2node = {nd["code"]: nd for nd in nodes}
    adj = {nd["code"]: [] for nd in nodes}
    for sg in segments:
        adj[sg["na"]].append(sg["id"])
        adj[sg["nb"]].append(sg["id"])
    seg_by_id = {sg["id"]: sg for sg in segments}

    def other_end(sg, code):
        return sg["nb"] if sg["na"] == code else sg["na"]

    unseen = set(sg["id"] for sg in segments)
    runs = []

    def walk(start_code, seg_id):
        path = []
        cur_code, cur_seg = start_code, seg_id
        while True:
            unseen.discard(cur_seg)
            path.append(cur_seg)
            nxt = other_end(seg_by_id[cur_seg], cur_code)
            if code2node[nxt]["degree"] != 2:
                return path, start_code, nxt
            cands = [s for s in adj[nxt] if s in unseen]
            if not cands:
                return path, start_code, nxt
            cur_seg, cur_code = cands[0], nxt

    for nd in sorted(nodes, key=lambda r: r["id"]):
        if nd["degree"] == 2:
            continue
        for sid in list(adj[nd["code"]]):
            if sid in unseen:
                path, na, nb = walk(nd["code"], sid)
                runs.append({"segs": path, "na": na, "nb": nb})
    while unseen:
        sid = next(iter(unseen))
        sg = seg_by_id[sid]
        path, na, nb = walk(sg["na"], sid)
        runs.append({"segs": path, "na": na, "nb": nb})
    for i, r in enumerate(runs, 1):
        r["run"] = f"走线{i}"
        r["length"] = sum(seg_by_id[s]["length"] for s in r["segs"])
    sc = 0
    for r in runs:
        for sid in r["segs"]:
            sc += 1
            seg_by_id[sid]["seg"] = f"SEG{sc:02d}"
            seg_by_id[sid]["run"] = r["run"]
    # 7) 输出组装
    segments_out = [{"seg": sg["seg"], "branch": branch_data[sg["branch"]]["key"],
                     "tube": branch_data[sg["branch"]]["occ"],
                     "na": sg["na"], "nb": sg["nb"],
                     "length": round(sg["length"], 1),
                     "dia": branch_data[sg["branch"]]["dia"],
                     "p0": [round(float(x), 1) for x in sg["p0"]],
                     "p1": [round(float(x), 1) for x in sg["p1"]],
                     "run": sg["run"]}
                    for sg in sorted(segments, key=lambda r: r["seg"])]
    nodes_out = [{"id": nd["id"], "code": nd["code"], "name_3d": nd["name_3d"],
                  "type": nd["type"],
                  "xyz": [round(float(x), 1) for x in nd["xyz"]],
                  "degree": nd["degree"],
                  "segments": sorted(set(seg_by_id[s]["seg"] for s, _e in nd["members"])),
                  "entity": nd["entity"], "nearest": nd["nearest"],
                  "secondary": [{"code": sc, "kind": sk, "tag": st}
                                for sc, sk, st in nd.get("secondary_codes", [])]}
                 for nd in nodes]
    entities = []
    for sg in sorted(segments, key=lambda r: r["seg"]):
        b = branch_data[sg["branch"]]
        entities.append({"code": sg["seg"], "category": "线段",
                         "name_3d": b["occ"],
                         "note": f"对应分支{b['key']}; 所属{sg['run']}"})
    for nd in nodes:
        if nd["type"] == "connector":
            extra = sorted(set(tag for k, tag in nd["entities"] if k == "connector"
                               and tag != nd["name_3d"]))
            note = f"位于节点{nd['code']}"
            if extra:
                note += "; 同位置实体:" + ",".join(extra)
            merged_from = sorted(k for k, v in relations.get("_tag_alias", {}).items()
                                 if v == nd["name_3d"])
            if merged_from:
                note += "; 共位合并:" + ",".join(merged_from)
            entities.append({"code": nd["code"], "category": "连接器",
                             "name_3d": nd["name_3d"], "note": note})
        elif nd["type"] == "clamp":
            extra = sorted(set(tag for k, tag in nd["entities"] if k == "clamp"
                               and tag != nd["name_3d"]))
            note = f"位于节点{nd['code']}"
            if extra:
                note += "; 同位置实体:" + ",".join(extra)
            merged_from = sorted(k for k, v in relations.get("_tag_alias", {}).items()
                                 if v == nd["name_3d"])
            if merged_from:
                note += "; 共位合并:" + ",".join(merged_from)
            entities.append({"code": nd["code"], "category": "固定卡扣",
                             "name_3d": nd["name_3d"], "note": note})
        elif nd["type"] == "tie":
            extra = sorted(set(tag for k, tag in nd["entities"] if k == "tie"
                               and tag != nd["name_3d"]))
            note = f"位于节点{nd['code']}"
            if extra:
                note += "; 同位置实体:" + ",".join(extra)
            entities.append({"code": nd["code"], "category": "扎带",
                             "name_3d": nd["name_3d"], "note": note})
        # 同点次要设备: 单独编码, 非拓扑节点
        for scode, sk, stag in nd.get("secondary_codes", []):
            cat = {"connector": "连接器", "clamp": "固定卡扣"}.get(sk, "扎带")
            merged_from = sorted(k for k, v in relations.get("_tag_alias", {}).items()
                                 if v == stag)
            note = f"与{nd['code']}同位置(次要设备, 非拓扑节点)"
            if merged_from:
                note += "; 共位合并:" + ",".join(merged_from)
            entities.append({"code": scode, "category": cat,
                             "name_3d": stag, "note": note})
    for nd in nodes:
        if nd["type"] in ("fork", "hanging", "joint"):
            entities.append({"code": nd["code"],
                             "category": "分支点" if nd["type"] == "fork" else "节点",
                             "name_3d": "",
                             "note": "自动命名; 相连线段" + ",".join(
                                 sorted(set(seg_by_id[s]["seg"] for s, _e in nd["members"])))})
    runs_out = [{"run": r["run"], "segs": [seg_by_id[s]["seg"] for s in r["segs"]],
                 "length": round(r["length"], 1), "na": r["na"], "nb": r["nb"]}
                for r in runs]
    branch_info = []
    for bi, b in enumerate(branch_data):
        sids = branch_seg[bi]
        branch_info.append({"key": b["key"],
                            "segs": [seg_by_id[s]["seg"] for s in sids],
                            "node0": seg_by_id[sids[0]]["na"],
                            "node1": seg_by_id[sids[-1]]["nb"]})
    return segments_out, nodes_out, entities, runs_out, branch_info


def _cluster_tubes(kept, log):
    """管段聚类成"分支": 端点相连且相切的拼成一条链;
    搭接式接头拆分(对标 VBA"端点落在他曲线身上=分支点"):
    链内两管段端头相接, 若接头落在第三条管体侧壁接触范围内,
    则此处为双分支搭接点, 不得合并为一根, 拆成独立分支.
    kept: [{spine, spine_len, solid:{vol,tag}, ...}]. 返回 clusters."""
    segs = [{"pts": it["spine"], "len": it["spine_len"], "it": it} for it in kept]
    pre_chains = list(chain_all(segs, tol=2.0, max_angle_deg=30.0))

    # 搭接式接头拆分(对标 VBA"端点落在他曲线身上=分支点"):
    # 链内两管段端头相接, 若接头落在第三条管体侧壁接触范围内,
    # 则此处为双分支搭接点, 不得合并为一根, 拆成独立分支
    chain_r = []
    for ch in pre_chains:
        rs = [_item_radius(m["it"]) for m in ch["members"]]
        chain_r.append(max(rs) if rs else 5.0)
    forbid = set()
    for ci, ch in enumerate(pre_chains):
        mems = ch["members"]
        if len(mems) < 2:
            continue
        for ai in range(len(mems)):
            for bi in range(ai+1, len(mems)):
                J = _spine_joint(mems[ai]["it"]["spine"], mems[bi]["it"]["spine"])
                if J is None:
                    continue
                for cj, ch2 in enumerate(pre_chains):
                    if cj == ci:
                        continue
                    pts2 = [tuple(p) for p in ch2["pts"]]
                    s, _q, dev = _poly_project(J, pts2)
                    if (dev <= chain_r[ci] + chain_r[cj] + CONTACT_TOL
                            and STATION_END_TOL < s < ch2["len"] - STATION_END_TOL):
                        log(f"搭接式接头(拆分不合并): {mems[ai]['it']['solid']['tag']} + "
                            f"{mems[bi]['it']['solid']['tag']} 接头落在管体站位{s:.1f}mm")
                        forbid.add((id(mems[ai]), id(mems[bi])))
                        break
    clusters = []
    for ch in pre_chains:
        mems = ch["members"]
        par = list(range(len(mems)))
        def _uf(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        for ai in range(len(mems)):
            for bi in range(ai+1, len(mems)):
                if (id(mems[ai]), id(mems[bi])) in forbid:
                    continue
                if _spine_joint(mems[ai]["it"]["spine"],
                                mems[bi]["it"]["spine"]) is None:
                    continue
                ra, rb = _uf(ai), _uf(bi)
                if ra != rb:
                    par[rb] = ra
        groups = defaultdict(list)
        for i, m in enumerate(mems):
            groups[_uf(i)].append(m)
        for g in groups.values():
            gs = [{"pts": m["it"]["spine"], "len": m["it"]["spine_len"],
                   "it": m["it"]} for m in g]
            for ch2 in chain_all(gs, tol=2.0, max_angle_deg=30.0):
                clusters.append({"items": [m["it"] for m in ch2["members"]],
                                 "spine": ch2["pts"], "spine_len": ch2["len"]})
    return clusters


def analyze(step_path, tol=3.0, out_dir=None, progress=None, force_reverse=False):
    def log(msg):
        if progress: progress(msg)

    doc = TDocStd_Document(TCollection_ExtendedString("XmlOcaf"))
    reader = STEPCAFControl_Reader()
    if reader.ReadFile(step_path) is None:
        raise RuntimeError(f"无法读取 STEP 文件: {step_path}")
    reader.Transfer(doc)
    st = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())

    tops = Sequence_TDF_Label(); st.GetShapes(tops)
    cands = []
    for i in range(1, tops.Length() + 1):
        lab = tops.Value(i)
        if S.IsAssembly_s(lab):
            cands.append((count_topo(S.GetShape_s(lab), TopAbs_FACE), lab))
    if not cands:
        raise RuntimeError("STEP 中没有找到装配结构")
    cands.sort(key=lambda t: -t[0])
    root, root_name = cands[0][1], label_name(cands[0][1])
    log(f"根装配: {root_name}")

    leaves = []
    def walk(lab, path, parent_trsf, proto_ctx):
        ref = TDF_Label()
        is_comp = S.IsComponent_s(lab)
        if is_comp:
            if not S.GetReferredShape_s(lab, ref):
                ref = lab
        else:
            ref = lab
        ref_name = label_name(ref)
        if is_real_name(ref_name):
            proto_ctx = ref_name
        if S.IsAssembly_s(ref):
            cs = Sequence_TDF_Label(); S.GetComponents_s(ref, cs)
            for i in range(1, cs.Length() + 1):
                c = cs.Value(i)
                loc = comp_location(c)
                t = loc.Transformation() if not loc.IsIdentity() else gp_Trsf()
                g = gp_Trsf(); g.Multiply(parent_trsf); g.Multiply(t)
                walk(c, path + [label_name(c)], g, proto_ctx)
        else:
            leaves.append({"occ": "/".join([p for p in path if is_real_name(p)]),
                           "proto": proto_ctx,
                           "shape": S.GetShape_s(ref),
                           "trsf": parent_trsf})

    walk(root, [], gp_Trsf(), root_name)
    log(f"装配实例数: {len(leaves)}")

    # 按产品原型收集: 线框 / 实体
    protos = defaultdict(lambda: {"wire": [], "wire_len": 0.0, "solids": [], "occ": ""})
    others = []
    for lf in leaves:
        shp = lf["shape"]
        nf = count_topo(shp, TopAbs_FACE)
        ne = count_topo(shp, TopAbs_EDGE)
        pd = protos[lf["proto"]]
        pd["occ"] = lf["occ"]
        if nf == 0 and ne > 0:
            curves, total = [], 0.0
            ex = TopExp_Explorer(shp, TopAbs_EDGE)
            while ex.More():
                e = TopoDS_Edge_s(ex.Current())
                L = edge_length(e); total += L
                curves.append({"len": L, "pts": xform_pts(sample_edge(e), lf["trsf"])})
                ex.Next()
            pd["wire"].append({"curves": curves, "occ": lf["occ"]})
            pd["wire_len"] += total
        elif nf > 0:
            # 扁平 STEP 常把多个实体做成一个复合体: 拆成单个实体逐个判定,
            # 否则整团实体会被 is_tube_solid 误判为非管状
            sols = []
            exs = TopExp_Explorer(shp, TopAbs_SOLID)
            while exs.More():
                sols.append(exs.Current())
                exs.Next()
            multi = len(sols) > 1
            base_tag = lf["occ"] or lf["proto"]
            for si, sol in enumerate(sols if multi else [shp]):
                sub = sol if multi else shp
                sub_nf = count_topo(sub, TopAbs_FACE)
                vol, bc = solid_volume_bbox(sub, lf["trsf"])
                tube = is_tube_solid(sub)
                tag = f"{base_tag}#S{si}" if multi else lf["occ"]
                pd["solids"].append({"occ": lf["occ"], "tag": tag, "proto": lf["proto"],
                                     "faces": sub_nf, "vol": vol, "bbox_c": bc,
                                     "shape": sub, "trsf": lf["trsf"], "tube": tube})
        else:
            others.append({"occ": lf["occ"], "proto": lf["proto"]})

    # 结构化分支判定(不依赖命名)
    # 1) 全部管状实体逐个反推中线
    tube_items = []
    for proto, pd in protos.items():
        for s in pd["solids"]:
            if not s["tube"]:
                continue
            r = reverse_spine(s["tube"], s["trsf"])
            if r:
                tube_items.append({"proto": proto, "solid": s,
                                   "spine": r[0], "spine_len": r[1]})
            else:
                log(f"警告: {s['tag']} 反推中心线失败, 跳过")
    # 2) 管状实体聚类成"分支": 同一产品内, 中线重叠的(半壳/分面)去重, 端点相连且
    #    相切的(同一扫掠体的分段)拼成一条; 分叉节点处有夹角的管子不合并

    clusters = []
    for proto, pd in protos.items():
        items = [it for it in tube_items if it["proto"] == proto]
        if not items:
            continue
        kept = []
        for it in sorted(items, key=lambda x: -x["spine_len"]):
            if any(spines_overlap(it["spine"], k["spine"], tol=2.0, frac=0.85)
                   for k in kept):
                continue
            kept.append(it)
        for cl in _cluster_tubes(kept, log):
            clusters.append(cl)
    # 3) 线框按产品分组成链(保留全部链, 逐条判定)
    wire_chains = []
    for proto, pd in protos.items():
        if pd["wire_len"] < WIRE_MIN_LEN:
            continue
        curves = [c for w in pd["wire"] for c in w["curves"]]
        for ch in chain_all(curves):
            if ch["len"] >= WIRE_MIN_LEN:
                wire_chains.append({"proto": proto, "pts": ch["pts"],
                                    "len": ch["len"], "n_seg": ch["n_seg"]})
    # 4) 每个聚类 -> 一个分支. 中线优先用空间匹配的线框(真中心线), 否则用实体反推;
    #    对不上任何分支的碎线框(如管表面的装饰线)直接忽略, 不再误判成分支
    branch_data = []
    used_solid_ids = set()
    for cl in clusters:
        items = cl["items"]
        protos_cl = set(it["proto"] for it in items)
        pts_rev, L_rev = cl["spine"], cl["spine_len"]
        vols = [it["solid"]["vol"] for it in items]
        vol = float(sum(vols))
        dia_est = 2*math.sqrt(vol/(math.pi*L_rev)) if L_rev > 0 else 0.0
        wire_hit = None
        if not force_reverse:
            for wc in wire_chains:
                if wc["proto"] not in protos_cl:
                    continue
                if not (0.5*L_rev <= wc["len"] <= 1.6*L_rev):
                    continue
                d = polyline_dist(wc["pts"][::max(1, len(wc["pts"])//60)], pts_rev)
                if np.mean(d) < max(2.0, 0.25*dia_est):
                    wire_hit = wc
                    break
        if wire_hit:
            pts, L, src, n_seg = wire_hit["pts"], wire_hit["len"], "wireframe", wire_hit["n_seg"]
        else:
            pts, L, src, n_seg = pts_rev, L_rev, "reversed", 0
        dia = round(2*math.sqrt(vol/(math.pi*L)), 2) if vol and L > 0 else None
        main_proto = max(protos_cl, key=lambda p: sum(1 for it in items if it["proto"] == p))
        occ = ";".join(sorted(set(it["solid"]["tag"] for it in items)))
        branch_data.append({"proto": main_proto, "occ": occ, "n_seg": n_seg,
                            "length": L, "vol": vol, "dia": dia,
                            "p0": pts[0], "p1": pts[-1], "pts": pts, "src": src,
                            "solids": [{"shape": it["solid"]["shape"], "trsf": it["solid"]["trsf"],
                                        "tag": it["solid"]["tag"], "bbox_c": it["solid"]["bbox_c"]}
                                       for it in items]})
        for it in items:
            used_solid_ids.add(id(it["solid"]))
        log(f"分支 {main_proto} [{occ}]: 长 {L:.1f}mm, 中心线来源={'线框' if src == 'wireframe' else '实体反推'}")
    # 5) 无实体、只有长线框的产品 -> 分支(旧行为保留)
    for proto, pd in protos.items():
        if pd["solids"] or pd["wire_len"] < WIRE_MIN_LEN:
            continue
        curves = [c for w in pd["wire"] for c in w["curves"]]
        pts, L = chain_curves(curves)
        if L < WIRE_MIN_LEN:
            continue
        branch_data.append({"proto": proto, "occ": pd["occ"], "n_seg": len(curves),
                            "length": L, "vol": None, "dia": None,
                            "p0": pts[0], "p1": pts[-1], "pts": pts, "src": "wireframe",
                            "solids": []})
        log(f"分支 {proto}: 长 {L:.1f}mm, 中心线来源=线框(纯线框)")
    # 6) 连接器候选: 未被分支消耗的实体(扁平文件里与分支同产品的连接器也能被收录)
    conn_candidates = [s for pd in protos.values() for s in pd["solids"]
                       if id(s) not in used_solid_ids]

    if not branch_data:
        raise RuntimeError("未识别到线束分支(无可用线框, 且无管状实体可反推)")

    # 分支唯一键: 同名分支(扁平文件)用实体编号 #S 区分; 名称唯一时保持原产品名
    _cnt = Counter(b["proto"] for b in branch_data)
    for b in branch_data:
        if _cnt[b["proto"]] > 1:
            tags = sorted(set(t.split("#S")[-1] for t in b["occ"].split(";") if "#S" in t),
                          key=lambda x: int(x) if x.isdigit() else x)
            b["key"] = b["proto"] + ("#S" + "+".join(tags) if tags else "")
        else:
            b["key"] = b["proto"]

    # ---- 接触式拓扑分析(2026-09-27 方案) ----
    try:
        relations = compute_relations(branch_data, conn_candidates, log)
    except Exception as e:
        log(f"接触分析失败({e}), 仅用端点邻近回退")
        relations = {"taps": [], "end_ends": [], "terminals": [], "tie_stations": []}
    # 端点邻近兜底(旧逻辑兼容): 端点距离<=容差即视为相连
    seen_ee = set()
    for (i, ei, j, ej) in relations["end_ends"]:
        seen_ee.add((i, ei, j, ej)); seen_ee.add((j, ej, i, ei))
    for i in range(len(branch_data)):
        for j in range(i+1, len(branch_data)):
            for ei, ej in ((0, 0), (0, 1), (1, 0), (1, 1)):
                if (i, ei, j, ej) in seen_ee:
                    continue
                pi = branch_data[i]["p0" if ei == 0 else "p1"]
                pj = branch_data[j]["p0" if ej == 0 else "p1"]
                if _d3(tuple(pi), tuple(pj)) <= tol:
                    relations["end_ends"].append((i, ei, j, ej))
    segments, nodes, entities, runs, branch_info = build_topology_from_relations(
        branch_data, relations, tol, conn_candidates, log)
    log(f"拓扑线段: {len(segments)}, 拓扑节点: {len(nodes)}, 走线: {len(runs)}")

    result = {
        "file": os.path.basename(step_path), "node_tol": tol,
        "branches": [], "segments": segments, "nodes": nodes,
        "entities": entities, "runs": runs, "connectors": [],
        "ignored_contacts": [{"a": x["a"], "b": x["b"],
                             "dist": round(x["dist"], 2),
                             "note": "侧-侧接触(两端均无管端进入对方实体), 未自动拆分, 需人工核对"}
                            for x in relations["side_sides"]],
    }
    for bi, b in enumerate(branch_data):
        pts = b["pts"]; step = max(1, len(pts)//120)
        result["branches"].append({
            "key": b["key"], "proto": b["proto"], "occ": b["occ"], "n_seg": b["n_seg"],
            "length": round(b["length"], 1), "dia": b["dia"], "src": b["src"],
            "p0": [round(float(x), 1) for x in b["p0"]],
            "p1": [round(float(x), 1) for x in b["p1"]],
            "polyline": [[round(float(x), 1) for x in p] for p in pts[::step]],
            "segs": branch_info[bi]["segs"],
            "node0": branch_info[bi]["node0"], "node1": branch_info[bi]["node1"]})
    result["connectors"] = [{"occ": s["tag"], "proto": s["proto"],
                             "center": [round(float(x), 1) for x in s["bbox_c"][:3]],
                             "bbox": [round(float(x), 1) for x in s["bbox_c"][3:]]}
                            for s in conn_candidates]

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "topology.json"), "w") as f:
            json.dump(result, f, indent=1)
    return result


if __name__ == "__main__":
    args = sys.argv[1:]
    force = "--force-reverse" in args
    args = [a for a in args if not a.startswith("--")]
    step = args[0] if len(args) > 0 else "M1E-DRD.stp"
    tol = float(args[1]) if len(args) > 1 else 3.0
    r = analyze(step, tol, out_dir=".", progress=print, force_reverse=force)
    print(f"\n分支: {len(r['branches'])}, 节点: {len(r['nodes'])}")
    for b in r["branches"]:
        print(f'{b["proto"]:24s} len={b["length"]:8.1f}mm dia~{b["dia"]}mm src={b["src"]}')
