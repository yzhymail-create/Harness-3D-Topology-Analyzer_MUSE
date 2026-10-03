#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""膨胀包围盒 + sweep-and-prune 接触检测优化

优化思路:
1. 膨胀AABB: 每个包围盒向外扩展 CONTACT_TOL (2mm)
2. sweep-and-prune: 对三个轴分别排序, 找出重叠的对
3. 只对这些候选对调用 BRepExtrema

复杂度: O(n log n) vs 原来的 O(n²)
"""

import numpy as np
from typing import List, Tuple, Dict, Any


def inflated_aabb(box, inflate: float) -> Tuple[float, float, float, float, float, float]:
    """膨胀AABB: 向外扩展 inflate mm
    
    Args:
        box: OCP.Bnd.Bnd_Box
        inflate: 膨胀距离 (mm)
    
    Returns:
        (xmin, ymin, zmin, xmax, ymax, zmax)
    """
    xmin = box.GetXMin() - inflate
    ymin = box.GetYMin() - inflate
    zmin = box.GetZMin() - inflate
    xmax = box.GetXMax() + inflate
    ymax = box.GetYMax() + inflate
    zmax = box.GetZMax() + inflate
    return (xmin, ymin, zmin, xmax, ymax, zmax)


def aabb_overlap(a: Tuple, b: Tuple) -> bool:
    """判断两个AABB是否重叠"""
    return not (a[3] < b[0] or b[3] < a[0] or
                a[4] < b[1] or b[4] < a[1] or
                a[5] < b[2] or b[5] < a[2])


def sweep_and_prune(items: List[Dict], key_fn, contact_tol: float) -> List[Tuple[int, int]]:
    """Sweep-and-prune 算法找出所有可能接触的实体对
    
    Args:
        items: 实体列表, 每个实体有 'box' 属性
        key_fn: 获取实体的包围盒的函数
        contact_tol: 接触容差 (mm)
    
    Returns:
        可能接触的实体对索引列表 [(i, j), ...]
    """
    n = len(items)
    if n < 2:
        return []
    
    # 提取膨胀后的AABB
    aabbs = [inflated_aabb(key_fn(item), contact_tol) for item in items]
    
    # 对X轴排序, 生成事件列表
    events = []
    for i, aabb in enumerate(aabbs):
        events.append((aabb[0], 'start', i))  # 左边界
        events.append((aabb[3], 'end', i))    # 右边界
    
    events.sort(key=lambda e: (e[0], 0 if e[1] == 'start' else 1))
    
    # Sweep: 维护活跃集合
    active = set()
    candidates = set()
    
    for pos, typ, idx in events:
        if typ == 'start':
            # 检查与所有活跃实体的Y/Z重叠
            for other in active:
                if aabb_overlap(aabbs[idx], aabbs[other]):
                    # 确保 i < j
                    pair = (min(idx, other), max(idx, other))
                    candidates.add(pair)
            active.add(idx)
        else:  # 'end'
            active.remove(idx)
    
    return list(candidates)


def test_sweep_and_prune():
    """测试 sweep-and-prune 算法"""
    from OCP.Bnd import Bnd_Box
    
    # 创建测试数据: 3个包围盒
    boxes = []
    
    # Box 0: (0,0,0) - (10,10,10)
    b0 = Bnd_Box()
    b0.Update(0, 0, 0, 10, 10, 10)
    boxes.append({'id': 0, 'box': b0})
    
    # Box 1: (8,8,8) - (20,20,20) - 与Box 0重叠
    b1 = Bnd_Box()
    b1.Update(8, 8, 8, 20, 20, 20)
    boxes.append({'id': 1, 'box': b1})
    
    # Box 2: (50,50,50) - (60,60,60) - 不与任何盒子重叠
    b2 = Bnd_Box()
    b2.Update(50, 50, 50, 60, 60, 60)
    boxes.append({'id': 2, 'box': b2})
    
    # 测试: contact_tol = 2mm
    candidates = sweep_and_prune(boxes, lambda x: x['box'], contact_tol=2.0)
    
    print(f"候选对: {candidates}")
    assert len(candidates) == 1
    assert candidates[0] == (0, 1)
    print("✓ 测试通过")


if __name__ == "__main__":
    test_sweep_and_prune()
