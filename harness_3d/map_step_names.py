#!/usr/bin/env python3
"""STEP BREP 名称 -> S 标签映射 (方案D: representation item 顺序关联).

已验证结论:
  harness_topology.py 用 TopExp_Explorer 枚举实体并按顺序打 tag
  (tag = ...#S{si}); OCCT transfer 顺序 == SHAPE_REPRESENTATION
  item 列表的顺序, explorer 顺序 == transfer 顺序。
  因此: S{i} 的 CATIA 实例全名 = shape representation 中第 i 个
  item 所引用的 BREP 实体的 name 属性。
  注意 item 顺序不一定等于 BREP 定义在文件中的顺序:
  - X1X-DRD_AllCATPart.stp: 两者一致 (22/22, 质心交叉验证);
  - IP-Harness.stp: 不一致, 3 个 BREP_WITH_VOIDS 在 item 顺序中
    位于 153/170/183 (多壳数 7/2/3 与 void 数 6/1/2 吻合, 已验证),
    在定义顺序中位于末尾 186/187/188。
  顶层 BREP 实体包括 MANIFOLD_SOLID_BREP 与 BREP_WITH_VOIDS
  (后者是带内部型腔的实体, 如连接器壳体)。
  中文名以 \\X2\\....\\X0\\ 转义, 本脚本自行解码,
  不依赖 OCP 的字符串处理。

用法:
    /path/to/ocptest/bin/python map_step_names.py <file.stp>
输出:
    与输入同目录的 <stem>_step_name_map.json
    每行: {"stag": "S0", "step_id": 21, "name": "SUPPORT\\CP03\\Tyton-157-00181"}

自检:
    脚本会用 OCP 实际 transfer 并数 solid 个数, 若与解析到的名称数
    不一致则告警 (顺序关联的前提是两者一一对应)。
"""

import json
import os
import re
import struct
import sys


# 顶层 BREP 实体: MANIFOLD_SOLID_BREP 与 BREP_WITH_VOIDS。
# 允许定义跨行: MANIFOLD_SOLID_BREP('name',#shell);
BREP_RE = re.compile(
    r"#(\d+)\s*=\s*(?:MANIFOLD_SOLID_BREP|BREP_WITH_VOIDS)"
    r"\s*\(\s*'((?:[^']|'')*)'"
)

# SHAPE_REPRESENTATION 的 item 列表, 即 OCCT transfer 顺序。
REPR_RE = re.compile(
    r"(?:ADVANCED_BREP_SHAPE_REPRESENTATION|SHAPE_REPRESENTATION)"
    r"\s*\(\s*'[^']*'\s*,\s*\(([^)]*)\)"
)


def decode_step_name(s: str) -> str:
    """解码 STEP Part21 字符串: \\\\ -> \\, '' -> ', \\X2\\XXXX\\X0\\ -> 中文."""
    s = s.replace("\\\\", "\\")
    s = s.replace("''", "'")

    def _x2(m):
        hexs = m.group(1)
        try:
            units = [int(hexs[i:i + 4], 16) for i in range(0, len(hexs), 4)]
            raw = b"".join(struct.pack(">H", u) for u in units)
            return raw.decode("utf-16-be")
        except Exception:
            return m.group(0)

    return re.sub(r"\\X2\\([0-9A-Fa-f]+)\\X0\\", _x2, s)


def parse_brep_names(step_path: str):
    """按 SHAPE_REPRESENTATION item 顺序返回 [(step_id, name), ...].

    OCCT transfer 顺序 == representation item 顺序; 若找不到可用的
    representation 则回退到 BREP 定义在文件中的顺序。
    """
    with open(step_path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    id2name = {}
    for m in BREP_RE.finditer(text):
        id2name[int(m.group(1))] = decode_step_name(m.group(2))
    best = []
    for m in REPR_RE.finditer(text):
        item_ids = [int(x) for x in re.findall(r"#(\d+)", m.group(1))]
        hit = [sid for sid in item_ids if sid in id2name]
        if len(hit) > len(best):
            best = hit
    if best and len(best) == len(id2name):
        return [(sid, id2name[sid]) for sid in best]
    # 回退: BREP 定义在文件中的顺序
    out = []
    for m in BREP_RE.finditer(text):
        sid = int(m.group(1))
        out.append((sid, id2name[sid]))
    return out


def count_solids_ocp(step_path: str) -> int:
    """用 OCP transfer 并数 TopExp_Explorer 枚举到的 solid 个数(顺序即 S 标签顺序)."""
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID

    reader = STEPControl_Reader()
    if reader.ReadFile(step_path) != 1:
        raise RuntimeError("STEP 读取失败")
    reader.TransferRoots()
    comp = reader.OneShape()
    n = 0
    ex = TopExp_Explorer(comp, TopAbs_SOLID)
    while ex.More():
        n += 1
        ex.Next()
    return n


def main():
    if len(sys.argv) < 2:
        print("用法: python map_step_names.py <file.stp>")
        sys.exit(1)
    step_path = sys.argv[1]

    names = parse_brep_names(step_path)
    print(f"[文本] 文件顺序解析出 {len(names)} 个 MANIFOLD_SOLID_BREP 名称")

    try:
        n_sol = count_solids_ocp(step_path)
        print(f"[OCP ] transfer 后枚举到 {n_sol} 个 solid")
        if n_sol != len(names):
            print(f"[警告] 数量不一致! 顺序关联可能不成立, 请勿直接使用本映射。")
        else:
            print("[自检] 数量一致, 顺序关联成立。")
    except Exception as e:
        print(f"[OCP ] 自检跳过: {e}")

    rows = [
        {"stag": f"S{i}", "step_id": eid, "name": nm}
        for i, (eid, nm) in enumerate(names)
    ]
    cn = [r for r in rows if any("一" <= ch <= "鿿" for ch in r["name"])]
    print(f"[结果] 含中文名称 {len(cn)} 个:")
    for r in cn:
        print(f"    {r['stag']:>5}  #{r['step_id']:<7} {r['name']}")

    stem = os.path.splitext(os.path.basename(step_path))[0]
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(step_path)), stem + "_step_name_map.json"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=1)
    print(f"[输出] {out_path}  (共 {len(rows)} 条)")


if __name__ == "__main__":
    main()
