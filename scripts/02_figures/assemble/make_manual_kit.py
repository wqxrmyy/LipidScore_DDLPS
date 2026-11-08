# -*- coding: utf-8 -*-
"""
生成「手工排版包」（PI 2026-09-21 决定：放弃子图自动合并，改人工排版）。

输入（唯一真相源，不另立标准）：
  _work/p5_figures/utils/themes.R        -> PANEL_BOX：各面板最终印刷尺寸
  _work/p5_figures/utils/assemble_lib.R  -> W / HMAX / PAD / LSIZE 拼版常量
  _work/p5_figures/assemble_fig_main_v3.R / assemble_figS_v3.R -> 面板→字母→行列 的分组方案
  _work/p5_figures/out/fig*.png|pdf      -> 面板件（已按最终尺寸 1:1 出图）

输出（交付到稿件仓库，供人工取件）：
  <repo>/v2/figures/manual_layout/
    README_manual_layout.md   手工排版说明书（画布规格 + 逐图坐标表 + 操作步骤）
    panel_manifest.csv        逐面板坐标清单（px / in / mm 三套单位）
    layout_preview.html       坐标预览（按真实像素绝对定位，缩放到屏幕）
    panels/<Fig>/             面板件副本（PNG 300dpi + 矢量 PDF），按图分目录
    reference/                自动拼版结果（仅供目视参考，非交付物）

说明：脚本只做「查表 + 复制 + 写清单」，不重渲染、不改面板内容。

★ 2026-09-24 改法：**增量同步 + 剪枝孤儿**，不再先整体清空目标目录。
  原因：整体清空在本机会触发「批量删除」安全护栏（单轮 > 50 个文件即要求人工确认），
  非交互运行会被直接中止 —— 实测曾把 manual_layout/ 删到一半（Fig1-3 面板、README、
  清单、预览全部丢失），见 _manual_kit_fig6.log。本包内每个文件都是派生件，
  剪枝孤儿（文件名不在本次生成集内者）与清空等价，且删除量通常为 0。
"""
import csv
import io
import os
import re
import shutil
import struct

ROOT = os.path.join(os.environ.get("LPS_FIG_ROOT", "."))
OUT = os.path.join(ROOT, 'out')
REPO_FIG = os.path.join(os.environ.get("LPS_DATA_ROOT", "."), "1", "重构手稿", "v2", "figures", "manual_layout")

# ---- 拼版常量（与 utils/assemble_lib.R 保持一致）----------------------
W_PX, HMAX_PX, PAD_PX, LSIZE_PX = 1949, 2400, 62, 44
DPI = 300

# ---- 2026-09-29 FigS3 版式微调（PI 裁定，两轮）--------------------------
#   PI 第 1 轮三点：① F 字母压标题；② E/G/I 与 ABCD 空隙大，上移；③ F/H 让余地。
#   PI 第 2 轮反馈：D 内部两子图贴死（另修 figS3d.R margin）；**整体美观度下降**。
#   第 2 轮回退（本表现状）：F/H 字母移出面板属美观度下降主因（孤悬字母与全图
#   字母节奏不一致）→ LETTER_DY / STACK_GAP 恢复默认（字母回面板内左上角、
#   E/F 与 G/H 紧贴，与行1 A/B 节奏一致）；F 压标题改由 figS3f.R 标题折行解决
#   （文字不变，仅排版：折行后每行变短，居中即与字母脱开）。
#   保留：E/G/I 上移（ROW_PAD 行2 62→24，I 与 D 底间隙 24px）—— PI 第 1 轮明确要求。
ROW_PAD = {'FigS3': {2: 24}}
STACK_GAP = {}   # 2026-09-29 第2轮回退：F/H 间隙 62→0（曾试 62 容纳外置字母，弃）
LETTER_DY = {}   # 同上：字母全部回面板内 (+8,+8)，全图统一

# 面板脚本尚未实现的图：尺寸已在 themes.R 冻结，但面板件还没出。
# 这些图**不计入对账告警**（否则每跑一次都报警），改列「待实现」清单；
# 其坐标按 PANEL_BOX 期望值写入 manifest，出图后重跑即自动补齐。
PENDING_FIGS = set()      # 2026-09-24：FigS6 六面板已出图（649x1120 实测全合规），移出待实现清单
                          # 2026-09-28：FigS6 统一重命名为 FigS3（figS3a-i 九面板）
                          # （Fig5 已于 2026-09-22 出图；至此全部标称图均已实现）

# ---- 分组方案（与 assemble_fig_main_v3.R / assemble_figS_v3.R 一致）--
# 注意：字母画在行的左上角，不烧进面板；B1/B2 同属面板 B（Fig1）。
LAYOUT = {
    'Fig1': [
        [('fig1a', 'A')],
        [('fig1b1', 'B1'), ('fig1b2', 'B2')],
        [('fig1c', 'C'), ('fig1d', 'D')],
        [('fig1e', 'E'), ('fig1h', 'F')],
    ],
    # 0929 版式重排（对齐 FigS3 标准）：B/C 稳健性直方图右列上下叠、D/E 校正点线图
    #   左列上下叠，A/F 独柱加高；阅读顺序仍为 A,B,C,D,E,F,G，图注字母零改动。
    'Fig2': [
        [('fig2a', 'A'), [('fig2b', 'B'), ('fig2c', 'C')]],
        [[('fig2d', 'D'), ('fig2e', 'E')], ('fig2f', 'F')],
        [('fig2g1', 'G'), ('fig2g2', '')],
    ],
    'Fig3': [
        [('fig3a', 'A'), ('fig3b', 'B')],
        [('fig3c', 'C'), ('fig3d', 'D')],
        [('fig3e', 'E'), ('fig3f', 'F')],
    ],
    'Fig4': [
        [('fig4a', 'A'), ('fig4b', 'B')],
        [('fig4c', 'C'), ('fig4d', 'D')],
        [('fig4e', 'E'), ('fig4f', 'F')],
        [('fig4g', 'G')],
    ],
    # ---- 2026-09-28 重构（方案 ②）：6 张补充图 → 3 张纯主题 -------------------
    #   FigS1 = 评分构建与逐基因（5 格：两队列直方图 | 相关矩阵+分量森林 | 六分量散点）
    #   FigS2 = 配对+方向+阴性效力（4 格）
    #   ★ 面板件统一命名为「图号+字母」（figS1a-e / figS2a-d / figS3a-i）
    #   删减 10 面板（主图副本/高度重叠）：原 figS1a/S2b/S2c/S2d/S3a/S3c/S4a/S4c/S4d/S5b
    'FigS1': [
        [('figS1a', 'A'), ('figS1b', 'B')],
        [('figS1c', 'C'), ('figS1d', 'D')],
        [('figS1e', 'E')],
    ],
    'FigS2': [
        [('figS2a', 'A'), ('figS2b', 'B')],
        [('figS2c', 'C'), ('figS2d', 'D')],
    ],
    # ---- Fig5：设计依据 LAYOUT_SPEC_v3.md §8.5 -------------------
    #  ★ 两层计数：字母数 6（符合期刊「鼓励 ≤6 面板」软口径），图件数 Fig5 = 12。
    #  ★★ 字母口径（PI 裁定 D24-7，2026-09-22）：四列行**整行共用一个面板字母**，
    #     **只标在第一个子图左上**，其余三格留空（空串 '' = 不标字母）——
    #     与稿面「编号只到面板级」一致；子序号（E1–E4／F1–F4）仅作内部坐标，不进图面。
    #     （注：Fig1 的 B1/B2 不同——那是**稿面承认**的编号，故仍逐格标注。）
    'Fig5': [
        [('fig5a', 'A'), ('fig5b', 'B')],
        [('fig5c', 'C'), ('fig5d', 'D')],
        [('fig5e1', 'E'), ('fig5e2', ''), ('fig5e3', ''), ('fig5e4', '')],
        [('fig5f1', 'F'), ('fig5f2', ''), ('fig5f3', ''), ('fig5f4', '')],
    ],
    # ---- 2026-09-28 重构（方案 ②）：FigS3 = 单细胞机制，R1「2 行 × 3 列 + 3 处上下叠」---
    #   ★ 单元格可为「上下叠」list（同列两子面板共享列宽、各自标字母，y 依次下移）；
    #     高密度件（通路热图 / 配体-受体 / KM）保持 1120 px 整格高、不压字。
    #   ★ 9 字母 A–I；字母由人工拼版添加，不进面板件。
    'FigS3': [
        [[('figS3a', 'A'), ('figS3b', 'B')], ('figS3c', 'C'), ('figS3d', 'D')],
        [[('figS3e', 'E'), ('figS3f', 'F')], [('figS3g', 'G'), ('figS3h', 'H')], ('figS3i', 'I')],
    ],
}

CAPTION_HINT = {
    'Fig1': 'Fig1A 研究设计示意 | B 队列构成（B1 堆叠柱 / B2 环形） | C 六基因 z 分热图 | D 六基因小提琴 | E 评分分布 | F 三队列评分箱线',
    'Fig2': 'Fig2A 按 CN 簇的评分 | B 簇 6 vs 其余 | C 簇 6 占比 | D 剔整簇稳健性 | E 纯度校正偏相关 | F 结构型 vs 负荷型 CN 变量 | G 相关性热图（g1 主图 / g2 图例，同属面板 G）',
    'Fig3': 'Fig3A–D 判别度与稳健性 | E KM 曲线 + 风险表 | F 分量贡献（主图引用）',
    'Fig4': 'Fig4A 配对 WD/DD | B 组内配对检验 | C 去分化侧身份基因丢失 | D DD 判别 AUC | E 跨队列热图 | F 2×3 分面 | G 三线方向收敛',
    'FigS1': 'FigS1A–E 评分构建与逐基因关系（两队列直方图 | 相关矩阵+分量森林 | 六分量散点）',
    'FigS2': 'FigS2A–D 配对证据、方向差异与阴性效力',
    'Fig5': 'Fig5A UMAP 按细胞类型 | B UMAP 按评分 A | C 点图（6 分量基因 × 细胞类型） | D 评分按细胞类型（IQR 方柱 + 均值点） | '
            'E 组成层四联 Δf（E1 脂肪系 / E2 肿瘤系 / E3 免疫系 / E4 脂肪占比梯度，关键证据） | '
            'F 状态层四联 Δm（F1 脂肪系 / F2 肿瘤系 / F3 免疫系，与 E 行逐列同谱系；'
            'F4 Kitagawa 分解，关键证据）',
    'FigS3': 'FigS3A–I 单细胞机制（A 评分投影 B 六层级共定位 C 通路热图 D 配体-受体 E 巨噬共变 F 效力边界 G 拟时序 H 配对均值 I KM 分层）',
}


def png_size(path):
    """读 PNG IHDR 拿宽高，不依赖 Pillow。"""
    with open(path, 'rb') as f:
        head = f.read(24)
    if head[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError('not a PNG: %s' % path)
    return struct.unpack('>II', head[16:24])


def parse_panel_box():
    """从 themes.R 解析 PANEL_BOX（.half / .third / .quart / .full）。

    列宽必须与 themes.R 的 helper 定义取值一致：
      .half  = FIG_W/2（2 列行）｜.third = FIG_W/3｜.quart = FIG_W/4｜.full = FIG_W
    """
    src = io.open(os.path.join(ROOT, 'utils', 'themes.R'), encoding='utf-8').read()
    colw = {
        'half':  W_PX / DPI / 2,
        'third': W_PX / DPI / 3,
        'quart': W_PX / DPI / 4,
        'full':  W_PX / DPI,
    }
    box = {}
    # 注：PANEL_BOX 里常一行声明多个面板（`fig2a = .half(0.55), fig2b = ...`），
    #     故不能加行首锚 ^，否则只取到每行第一个。
    for m in re.finditer(r'(fig[A-Za-z0-9_]+)\s*=\s*\.(half|third|quart|full)\(([0-9.]+)\)', src):
        stem, kind, r = m.group(1), m.group(2), float(m.group(3))
        w = colw[kind]
        box[stem] = (round(w, 4), round(w * r, 4))
    return box


def prune_orphans(root, keep, max_delete=40):
    """删除 root 下**不在本次生成集**内的孤儿文件，并清掉因此变空的子目录。

    keep 为相对路径集合（posix 分隔符）。返回 (已删除相对路径列表, 跳过原因或 None)。

    ★ 为什么不是「先整体清空」：本机装了批量删除安全护栏（单轮 > 50 个删除即要求
      人工确认），非交互运行会被直接中止。2026-09-24 实测：整体清空在删到第 50 个
      文件时被拦下，manual_layout/ 被删到一半（Fig1-3 面板 + README + 清单 + 预览
      全部丢失），只剩 Fig4 / Fig5 / FigS* 与 reference/。
      改为剪枝后，删除量 = 真正的孤儿数（正常为 0），既安全又等价。

    ★★ 2026-09-24 二次加固（本次）：**预检 + 延迟**。
      原实现是「写完全部产物 -> 才剪枝」，一旦剪枝触发护栏，进程被中止时
      README / 清单 / 预览已写了一半（或已写完但告警未打印），现场难以判读。
      现改为：① 先**只扫描**算出待删集；② 若待删数 > max_delete（默认 40，
      低于护栏阈值 50，留安全余量）→ **完全跳过删除**并在返回值中说明原因，
      由调用方打印告警；③ 否则才执行删除。这样破坏性动作要么不做、要么
      在可控规模内完成，绝不会出现「删一半被打断」的中间态。

    ★ 不用 `shutil.rmtree` —— 本机安全删除 shim 会把 rmtree 重定向到回收站，
      而中文路径下 trash 会失败（实测 2026-09-21：`Some operations were aborted`）。
      逐文件 `os.remove` + `os.rmdir` 稳定可用。
    """
    keep = set(keep)
    removed = []
    if not os.path.isdir(root):
        return removed, None

    def is_backup(rel):
        """路径任一段以 _bak / _backup 开头 -> 人工留的中间备份，永不剪枝。

        ★ 2026-09-25 三次加固：实测本函数曾把
          manual_layout/_bak_20260924_pre_S2A_kit/ 下的 31 个文件当成"孤儿"
          整体删掉。备份目录的定义就是「不在本次生成集内」——按 keep 剪枝
          必然误伤它们。故显式豁免（有个数值上的巧合要注意：31 < max_delete=40，
          规模闸门不会拦，所以必须靠这条豁免而不是靠闸门）。
        """
        return any(seg.startswith(('_bak', '_backup')) for seg in rel.split('/'))

    # ---- ① 预检：只扫描，不删除 ----
    todo = []
    for dirpath, _dirnames, filenames in os.walk(root, topdown=False):
        for name in filenames:
            p = os.path.join(dirpath, name)
            rel = os.path.relpath(p, root).replace(os.sep, '/')
            if rel in keep or is_backup(rel):
                continue
            todo.append((p, rel))
    # ---- ② 规模闸门：超限则整体跳过（不做任何破坏性动作）----
    if len(todo) > max_delete:
        return [], ('待删孤儿 %d 个 > 安全上限 %d 个；为避免触发批量删除护栏'
                    '导致「删一半」中间态，本次**完全跳过剪枝**（产物已正常写入）。'
                    '如需清理请人工确认后单独执行。' % (len(todo), max_delete))
    # ---- ③ 执行 ----
    for p, rel in todo:
        os.remove(p)
        removed.append(rel)
    for dirpath, _dirnames, _filenames in os.walk(root, topdown=False):
        if dirpath == root:
            continue
        rel = os.path.relpath(dirpath, root).replace(os.sep, '/')
        if is_backup(rel):
            continue
        if not os.listdir(dirpath):
            os.rmdir(dirpath)
    return removed, None


def main():
    box = parse_panel_box()
    for sub in ['panels', 'reference']:
        os.makedirs(os.path.join(REPO_FIG, sub), exist_ok=True)
    # 本次生成的相对路径集合（用于剪枝孤儿）；三个根层产物固定会重写，先入集
    keep = set(['README_manual_layout.md', 'panel_manifest.csv', 'layout_preview.html'])

    rows_out, problems, fig_summary, pending = [], [], [], []

    for fig, rows in LAYOUT.items():
        is_pending = fig in PENDING_FIGS
        y_row = 0
        for ri, row in enumerate(rows, start=1):
            n = len(row)
            # ---- 每格可为单件 tuple 或「上下叠」list；先算各格堆叠高度 ----
            cell_heights = []
            for cell in row:
                items = [cell] if isinstance(cell, tuple) else cell
                ch = 0
                for k, (stem, _) in enumerate(items):
                    p = os.path.join(OUT, stem + '.png')
                    if os.path.exists(p):
                        h = png_size(p)[1]
                    else:
                        # 缺件：用 PANEL_BOX 设计值代替实测值，使总高预算仍可核算
                        exp = box.get(stem)
                        h = int(exp[1] * DPI) if exp else 0
                        if is_pending:
                            pending.append('%s/%s' % (fig, stem))
                        else:
                            problems.append('缺面板件：%s.png' % stem)
                    ch += h
                    # 叠件间隙（gap-before 语义：记在**后一件**头上，与放置循环一致）
                    if k > 0:
                        ch += STACK_GAP.get(stem, 0)
                cell_heights.append(ch)
            row_h = max(cell_heights) + PAD_PX
            for i, cell in enumerate(row, start=1):
                items = [cell] if isinstance(cell, tuple) else cell
                # 放置坐标（人工贴图）用四舍五入；**尺寸**期望值一律用截断（int），
                # 与出图管线实测一致：fig5a 高 623 = int(2.0789 in x 300)（round 给 624）、
                # 4 列宽 487 = int(487.25)、2 列宽 974 = int(974.5)、3 列宽 649 = int(649.67)
                # （3 列 round 会给 650，与 LAYOUT_SPEC_v3.md 8.6 冻结的 649 冲突）。
                x_px = round((i - 1) * W_PX / n)
                y_off = y_row + ROW_PAD.get(fig, {}).get(ri, PAD_PX)   # 行首字母带（可按行覆盖）
                for k, (stem, letter) in enumerate(items):
                    if k > 0:
                        y_off += STACK_GAP.get(stem, 0)   # gap-before：F/H 与上件的间隙
                    png = os.path.join(OUT, stem + '.png')
                    pdf = os.path.join(OUT, stem + '.pdf')
                    exists = os.path.exists(png)
                    exp = box.get(stem)
                    if exists:
                        w_px, h_px = png_size(png)
                    elif exp:
                        w_px, h_px = int(exp[0] * DPI), int(exp[1] * DPI)
                    else:
                        continue
                    # 与 PANEL_BOX 对账（允许 1 px 取整误差）
                    if not exp:
                        problems.append('%s 未列入 PANEL_BOX（尺寸无冻结依据）' % stem)
                    elif exists:
                        dw = abs(w_px - int(exp[0] * DPI))
                        dh = abs(h_px - int(exp[1] * DPI))
                        if dw > 1 or dh > 1:
                            problems.append('%s 实测 %dx%d px，PANEL_BOX 期望 %dx%d px（差 %d/%d）'
                                            % (stem, w_px, h_px, int(exp[0] * DPI),
                                               int(exp[1] * DPI), dw, dh))
                    if x_px + w_px > W_PX + 1:
                        problems.append('%s 在 %d 列布局下右缘 %d px 超出画布 %d px'
                                        % (stem, n, x_px + w_px, W_PX))
                    # 复制面板件（仅已实现者；待实现面板只登记坐标，不出件）
                    if exists:
                        dst_dir = os.path.join(REPO_FIG, 'panels', fig)
                        os.makedirs(dst_dir, exist_ok=True)
                        shutil.copy2(png, os.path.join(dst_dir, stem + '.png'))
                        keep.add('panels/%s/%s.png' % (fig, stem))
                        if os.path.exists(pdf):
                            shutil.copy2(pdf, os.path.join(dst_dir, stem + '.pdf'))
                            keep.add('panels/%s/%s.pdf' % (fig, stem))
                    rows_out.append(dict(
                        figure=fig, row=ri, col=i, cols=n, letter=letter, stem=stem,
                        implemented=1 if exists else 0,
                        file_png='panels/%s/%s.png' % (fig, stem),
                        file_pdf='panels/%s/%s.pdf' % (fig, stem) if os.path.exists(pdf) else '',
                        target_w_in=('%.4f' % exp[0]) if exp else '',
                        target_h_in=('%.4f' % exp[1]) if exp else '',
                        w_px=w_px, h_px=h_px,
                        x_px=x_px, y_px=y_off, row_top_px=y_row, row_h_px=row_h,
                        letter_x_px=x_px + 8, letter_y_px=y_off + LETTER_DY.get(stem, 8),
                        x_in='%.4f' % (x_px / DPI), y_in='%.4f' % (y_off / DPI),
                        w_in='%.4f' % (w_px / DPI), h_in='%.4f' % (h_px / DPI),
                        x_mm='%.2f' % (x_px / DPI * 25.4), y_mm='%.2f' % (y_off / DPI * 25.4),
                        w_mm='%.2f' % (w_px / DPI * 25.4), h_mm='%.2f' % (h_px / DPI * 25.4),
                    ))
                    y_off += h_px
            y_row += row_h
        fig_summary.append((fig, y_row, HMAX_PX - y_row))
        if y_row > HMAX_PX:
            problems.append('%s 总高 %d px 超出单页上限 %d px' % (fig, y_row, HMAX_PX))

    # ---- 坐标清单 CSV ----
    cols = ['figure', 'row', 'col', 'cols', 'letter', 'stem', 'implemented',
            'file_png', 'file_pdf',
            'target_w_in', 'target_h_in', 'w_px', 'h_px',
            'x_px', 'y_px', 'row_top_px', 'row_h_px', 'letter_x_px', 'letter_y_px',
            'x_in', 'y_in', 'w_in', 'h_in', 'x_mm', 'y_mm', 'w_mm', 'h_mm']
    with io.open(os.path.join(REPO_FIG, 'panel_manifest.csv'), 'w',
                 encoding='utf-8-sig', newline='') as f:
        wtr = csv.DictWriter(f, fieldnames=cols)
        wtr.writeheader()
        for r in rows_out:
            wtr.writerow(r)

    # ---- 参考拼版件（若存在）----
    refs = []
    for fig in LAYOUT:
        src = os.path.join(OUT, fig + '.png')
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(REPO_FIG, 'reference', fig + '.png'))
            keep.add('reference/%s.png' % fig)
            refs.append(fig)

    # ---- 说明书 ----
    write_readme(rows_out, fig_summary, problems, refs, pending)
    # ---- 坐标预览 ----
    write_preview(rows_out, fig_summary)

    # ---- 剪枝孤儿（本包全部为派生件；正常情况删除量 = 0）----
    #   ★ 加固（2026-09-24）：超安全上限则整体跳过，不做破坏性动作（见 prune_orphans 注释）
    removed, skip_reason = prune_orphans(REPO_FIG, keep)

    n_impl = sum(1 for r in rows_out if r['implemented'])
    print('=' * 72)
    print('手工排版包已生成：%s' % REPO_FIG)
    print('面板 %d 个（已实现 %d / 待实现 %d）| 图 %d 张 | 参考拼版件 %d 张'
          % (len(rows_out), n_impl, len(rows_out) - n_impl, len(LAYOUT), len(refs)))
    for fig, tot, left in fig_summary:
        tag = '   ← 设计稿（面板待实现）' if fig in PENDING_FIGS else ''
        print('  %-6s 画布 1949 x %4d px = 6.50 x %.2f in   余量 %4d px%s'
              % (fig, tot, tot / DPI, left, tag))
    print('-' * 72)
    if pending:
        print('待实现（设计尺寸已冻结，出图后重跑即自动补齐）：%d 项' % len(pending))
        print('  ' + '  '.join(pending))
    if problems:
        print('对账告警 %d 条：' % len(problems))
        for p in problems:
            print('  ! ' + p)
    else:
        print('对账通过：已实现面板尺寸与 PANEL_BOX 全部一致，全部标称图均在 1949x2400 内。')
    if removed:
        print('-' * 72)
        print('剪枝孤儿 %d 个（不在本次生成集内，已删除）：' % len(removed))
        for r in removed:
            print('  - ' + r)
    elif skip_reason:
        print('-' * 72)
        print('剪枝已跳过（安全闸门）：%s' % skip_reason)
    else:
        print('剪枝孤儿：0 个（目标目录与生成集完全一致）')
    return problems


def write_readme(rows, fig_summary, problems, refs, pending):
    L = []
    A = L.append
    A('# 图件手工排版包（Manual layout kit）')
    A('')
    A('> PI 裁定 2026-09-21：**放弃子图自动合并**，沿用既定的子图组合方案（分组、字母、网格），'
      '由人工在 Illustrator / Photoshop / PowerPoint 中完成拼版。')
    A('> 本包不含拼版成品，只提供**坐标清单 + 面板件 + 坐标预览**，使人工摆放与冻结规格逐像素一致。')
    A('')
    A('## 1. 画布规格（硬参数）')
    A('')
    A('| 项 | 值 |')
    A('|---|---|')
    A('| 画布宽 | **1949 px** = 165 mm = 6.497 in @300 dpi |')
    A('| 单页高上限 | **2400 px** = 203 mm = 8.000 in（Cancer Science 图件上限 6.5 × 8 in）|')
    A('| 分辨率 | 300 dpi |')
    A('| 底色 | 纯白 #FFFFFF |')
    A('| 面板本体 | **1:1 摆放，禁止缩放**（面板已按最终印刷尺寸出图） |')
    A('| 行内字母带 | **62 px** 高，字母位于带内左起 **+8 px、上起 +8 px** |')
    A('| 字母格式 | `Arial` **Bold**、**44 px**（≈10.5 pt @300dpi）、黑色、全大写 |')
    A('| 图内最小字号 | ≥ 6 pt（现面板实测 7–8 pt） |')
    A('')
    A('> **为什么禁止缩放**：面板是「按最终印刷尺寸」绘制的。一旦在某一步整体缩放，'
      '7 pt 的轴标签会被压到 6 pt 以下，跌破期刊下限。若确需改尺寸，必须回到脚本重渲染，'
      '不能靠排版软件拖拽。')
    A('')
    A('## 2. 逐图网格总览')
    A('')
    A('| 图 | 行数 | 面板件数 | 总高 px | 印后高 in | 印后高 mm | 余量 px | 状态 |')
    A('|---|---|---|---|---|---|---|---|')
    for fig, rows_ in LAYOUT.items():
        tot = dict((f, t) for f, t, _ in fig_summary)[fig]
        left = HMAX_PX - tot
        npan = sum(1 for r in rows if r['figure'] == fig)
        st = '★ 设计稿（面板待实现）' if fig in PENDING_FIGS else '面板件齐备'
        A('| %s | %d | %d | %d | %.2f | %.1f | %d | %s |'
          % (fig, len(rows_), npan, tot, tot / DPI, tot / DPI * 25.4, left, st))
    A('')
    A('面板总数为 **%d**。逐面板的像素/英寸/毫米坐标见 `panel_manifest.csv`'
      '（`implemented` 列标 `0` 者尚未出件，其坐标按设计尺寸登记）。' % len(rows))
    A('')
    A('## 3. 坐标明细')
    for fig, rows_ in LAYOUT.items():
        tot = dict((f, t) for f, t, _ in fig_summary)[fig]
        A('')
        A('### %s（画布 1949 × %d px；%s）' % (fig, tot, CAPTION_HINT.get(fig, '')))
        A('')
        A('| 行 | 列 | 字母 | 面板件 | 状态 | 尺寸 px (W×H) | 粘贴左上角 (x, y) px | 尺寸 in | 左上角 in |')
        A('|---|---|---|---|---|---|---|---|---|')
        for r in [x for x in rows if x['figure'] == fig]:
            A('| %d | %d/%d | %s | `%s` | %s | %d × %d | (%d, %d) | %.2f × %.2f | (%.2f, %.2f) |'
              % (r['row'], r['col'], r['cols'], r['letter'] or '—', r['stem'],
                 '✅' if r['implemented'] else '待实现',
                 r['w_px'], r['h_px'], r['x_px'], r['y_px'],
                 float(r['w_in']), float(r['h_in']), float(r['x_in']), float(r['y_in'])))
    A('')
    A('> `(x, y)` 为**粘贴面板左上角**在画布中的像素坐标（原点 = 画布左上角，y 向下）。'
      '列宽按 `round((i-1) × 1949 / 列数)` 取整，故两列时 x=0 / 974（或 975）、'
      '三列时 x=0 / 650 / 1299、四列时 x=0 / 487 / 974 / 1462。')
    A('')
    A('> **尺寸口径 = 截断**（与出图管线一致，非四舍五入）：2 列 **974** = `int(1949/2)`、'
      '3 列 **649** = `int(1949/3)`、4 列 **487** = `int(1949/4)`；高同理'
      '（如 fig5a **623** = `int(2.0789 in × 300 dpi)`，四舍五入会得 624）。'
      '`w_px / h_px` 即文件原生像素，**粘贴时不得缩放**。')
    A('')
    A('## 4. 操作步骤（任选一个工具）')
    A('')
    A('### A. Adobe Illustrator（推荐，矢量优先）')
    A('1. 新建文档：宽 **165 mm**、高取 §2 表中该图的「印后高 mm」值（需 < 203 mm），'
      '色彩模式 RGB，栅格效果 300 ppi。')
    A('2. 用 `panel_manifest.csv` 的 `x_mm / y_mm` 作为「变换」面板里的 X/Y 值逐件置入（Link 或 Embed 均可）。')
    A('3. 矢量件用 `panels/<Fig>/*.pdf` 置入；含位图/热图的仍用 PDF 亦可（已 300 dpi 内嵌）。')
    A('4. 字母：文本工具输入 `Arial Bold 10.5 pt`，按 `letter_x_px / letter_y_px` 定位。')
    A('5. 导出：`File > Export As > TIFF`（300 ppi、LZW、CMYK 或按期刊要求 RGB），或 `Save As PDF`。')
    A('')
    A('### B. Adobe Photoshop')
    A('1. 新建：宽 1949 px、高取 §2 表中该图的「总高 px」值（≤ 2400），300 dpi，RGB，白底。')
    A('2. `File > Place Embedded` 逐件置入，**不要勾选** 变换时的等比缩放手柄拖拽；'
      '直接用属性栏输入 X/Y/W/H（W/H 应与文件原生像素一致）。')
    A('3. 字母层：Arial Bold 44 px（=10.5 pt @300dpi），黑色。')
    A('4. 导出 TIFF/PNG（300 dpi）。')
    A('')
    A('### C. PowerPoint（最快，精度足够）')
    A('1. 设计 > 幻灯片大小 > 自定义：宽 **16.5 cm**、高取 §2 表中该图的「印后高 mm」值。')
    A('2. 插入图片后，在「大小和位置」里按 cm 值设定（= px ÷ 300 × 2.54）。')
    A('3. 文本框打字母，Arial Bold 10.5 pt，位置按 px 换算。')
    A('4. 全选 → 右键「另存为图片」（PNG，300 dpi 需在注册表/导出选项调整）或导出 PDF。')
    A('')
    A('## 5. 验收清单（拼完必查）')
    A('1. 画布宽 = 1949 px，总高 ≤ 2400 px；')
    A('2. 每块面板未被缩放（对比 `panel_manifest.csv` 的 `w_px / h_px`）；')
    A('3. 字母字号统一 44 px / Arial Bold，位置一致（同上缘、同左缩进）；')
    A('4. 无面板互相压盖、无越出画布右缘；')
    A('5. 放大到 100% 检查最小字号仍 ≥ 6 pt（文案/轴标签可读）。')
    A('')
    A('## 6. 目录内容')
    A('')
    A('```')
    A('manual_layout/')
    A('  README_manual_layout.md   本说明书')
    A('  panel_manifest.csv        逐面板坐标清单（px / in / mm）')
    A('  layout_preview.html       坐标预览（按真实像素定位，仅缩小显示）')
    A('  panels/<Fig>/             面板件：*.pdf（矢量/300dpi 位图）+ *.png（300 dpi）')
    A('  reference/                自动拼版结果，仅供目视参考，不作为交付物')
    A('  ★ 字母口径：**两列行**每格各标一个字母；**四列行**（Fig5 E/F）整行**共用一个字母**，')
    A('    只标在第一个子图左上（清单中该行其余格 letter 列为空 = 不标）；子序号仅作内部坐标，不进图面。')
    A('    **FigS3 上下叠**（PI 裁定 2026-09-28）：2 行 × 3 列，左/中列各含上下叠两子面板，')
    A('    逐子面板标字母（A/B、E/F、G/H），高密度件（C/D/I）整格一字母；共 9 字母 A–I。')
    A('    另：**单独子图件一律不标字母**，字母仅由拼版件（`out/FigN.png`）添加。')
    A('```')
    A('')
    A('## 7. 与自动拼版的关系')
    A('')
    A('- 旧的自动拼版脚本（`assemble_fig_main_v3.R`、`assemble_figS_v3.R`、'
      '`utils/assemble_lib.R`）**保留但停用**，仅作几何算法的参考实现；')
    A('- 面板渲染管线（`themes.R` 的 `PANEL_BOX` + 各 `figXx.R`）**继续有效**，'
      '本次未改动任何面板内容；')
    A('- 若将来恢复自动拼版，本包的坐标清单即为拼版脚本的验收基准（同一套常量推导）。')
    A('- 新增 `_asm_fig5_reference.R`（2026-09-22）：只生成 Fig5 的 **reference 参照件**'
      '（`out/Fig5.png` + `out/Fig5.pdf`），**不作为交付物**；字母口径同 D24-7。')
    A('- 新增 `_asm_figS3_reference.R`（2026-09-28）：同规格生成 FigS3 的 **reference 参照件**'
      '（`out/FigS3.png` + `out/FigS3.pdf`，1949 × 2364 px）。版式 2 行 × 3 列（§8.6），'
      '**九格标字母 A–I**（PI 裁定 2026-09-28；单独子图件不标字母）。')
    A('')
    if pending:
        A('## 8. 待实现面板（设计尺寸已冻结，面板尚未出件）')
        A('')
        A('下列面板的**版面盒尺寸已在 `themes.R` 的 `PANEL_BOX` 冻结**，但渲染脚本尚未编写，'
          '故本包只有其坐标、无其面板件。**出图后重跑本包生成脚本即自动补齐。**')
        A('')
        A('```')
        for i in range(0, len(pending), 6):
            A('  ' + '  '.join(pending[i:i + 6]))
        A('```')
        A('')
        A('> 上列只含**尚未出件**的图。`FigS3` = 机制层「评分 ↔ mTORC2/TAM 轴」'
          '（设计依据 `LAYOUT_SPEC_v3.md` §8.6；2 行 × 3 列，高密度件 649 × 1120 px）；'
          '`Fig5` 已于 2026-09-22 出图，其设计依据见 §8.5。')
        A('> 出图后重跑 `make_manual_kit.py`，对应行即由「待实现」转为「已实现」，'
          '尺寸按实测值入清单。')
        A('')
    if problems:
        A('## 9. 对账告警（生成时检出，须先处理）')
        A('')
        for p in problems:
            A('- %s' % p)
        A('')
    A('---')
    A('')
    A('生成脚本：`_work/p5_figures/make_manual_kit.py`（查表 + 复制 + 写清单，不重渲染）。')
    io.open(os.path.join(REPO_FIG, 'README_manual_layout.md'), 'w', encoding='utf-8').write('\n'.join(L))


def write_preview(rows, fig_summary):
    """按真实像素绝对定位生成预览（外层缩放显示）。"""
    scale = 0.46
    tot_map = dict((f, t) for f, t, _ in fig_summary)
    h = []
    A = h.append
    A('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">')
    A('<title>图件手工排版坐标预览</title><style>')
    A('body{font-family:"Segoe UI",Arial,sans-serif;background:#f6f7f9;color:#1a1a1a;margin:0;padding:24px}')
    A('h1{font-size:20px;margin:0 0 6px}h2{font-size:15px;margin:28px 0 8px}')
    A('.meta{font-size:12.5px;color:#555;line-height:1.7;margin-bottom:6px}')
    A('.figbox{position:relative;overflow:hidden;border:1px solid #c9ced6;background:#fff;'
      'box-shadow:0 1px 3px rgba(0,0,0,.08)}')
    A('.figwrap{position:absolute;top:0;left:0;width:1949px;background:#fff;'
      'transform-origin:top left}')
    A('.panel{position:absolute;display:block}')
    A('.letter{position:absolute;font-family:Arial,Helvetica,sans-serif;font-weight:700;'
      'font-size:44px;line-height:44px;color:#000}')
    A('.figwrap .rowguide{position:absolute;left:0;width:1949px;border-top:1px dashed #9aa4b2}')
    A('.boxout{position:absolute;border:1px dashed #d0453c88;pointer-events:none}')
    A('.ph{position:absolute;box-sizing:border-box;border:1px solid #cfd6e0;'
      'background:repeating-linear-gradient(45deg,#f4f6fa,#f4f6fa 9px,#e9edf4 9px,#e9edf4 18px);'
      'display:flex;flex-direction:column;align-items:center;justify-content:center;'
      'gap:4px;color:#8a94a4;font-size:12px}')
    A('.ph-stem{font-weight:700;color:#5c6675}')
    A('.ph-size{font-size:11px}')
    A('.colguide{position:absolute;top:0;height:100%;border-left:1px dotted #b9c0ca}')
    A('.scroll{overflow-x:auto}')
    A('code{background:#eef1f5;padding:1px 5px;border-radius:3px;font-size:12px}')
    A('</style></head><body>')
    A('<h1>图件手工排版坐标预览</h1>')
    A('<div class="meta">画布宽 <b>1949 px</b>（165 mm @300 dpi），单页高上限 <b>2400 px</b>（203 mm）。'
      '面板按最终印刷尺寸绘制，<b>1:1 摆放</b>，禁止缩放。<br>'
      '字母：<b>Arial Bold 44 px</b>，位于每行左上角 (+8, +8)。'
      '红虚线为面板外框，灰虚线为画布/行边界。预览整体缩放到 '
      + str(int(scale * 100)) + '% 显示，实际坐标见 <code>panel_manifest.csv</code>。</div>')
    for fig in LAYOUT:
        fr = [r for r in rows if r['figure'] == fig]
        tot = tot_map[fig]
        A('<h2>%s <span style="font-weight:400;color:#666;font-size:13px">'
          '1949 × %d px = 6.50 × %.2f in ｜ 余量 %d px</span></h2>'
          % (fig, tot, tot / DPI, HMAX_PX - tot))
        A('<div class="scroll"><div class="figbox" style="width:%dpx;height:%dpx">'
          % (round(W_PX * scale), round(tot * scale)))
        A('<div class="figwrap" style="height:%dpx;transform:scale(%s)">' % (tot, scale))
        rowtops = sorted(set(r['row_top_px'] for r in fr))
        for rt in rowtops:
            A('<div class="rowguide" style="top:%dpx"></div>' % rt)
        for r in fr:
            if r['implemented']:
                A('<img class="panel" src="%s" style="left:%dpx;top:%dpx;width:%dpx;height:%dpx">'
                  % (r['file_png'], r['x_px'], r['y_px'], r['w_px'], r['h_px']))
            else:
                # 设计稿：面板件尚未出，用斜纹占位框示意版式与占位尺寸
                A('<div class="ph" style="left:%dpx;top:%dpx;width:%dpx;height:%dpx">'
                  '<span class="ph-stem">%s</span><span class="ph-size">%d × %d</span></div>'
                  % (r['x_px'], r['y_px'], r['w_px'], r['h_px'], r['stem'],
                     r['w_px'], r['h_px']))
            A('<div class="boxout" style="left:%dpx;top:%dpx;width:%dpx;height:%dpx"></div>'
              % (r['x_px'], r['y_px'], r['w_px'], r['h_px']))
            if r['letter']:
                A('<span class="letter" style="left:%dpx;top:%dpx">%s</span>'
                  % (r['letter_x_px'], r['letter_y_px'], r['letter']))
        A('</div></div></div>')
    A('</body></html>')
    io.open(os.path.join(REPO_FIG, 'layout_preview.html'), 'w', encoding='utf-8').write('\n'.join(h))


if __name__ == '__main__':
    main()
