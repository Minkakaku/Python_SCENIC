#!/usr/bin/env python
"""pySCENIC 独立分析入口：h5ad → GRNBoost2 → motif pruning → AUCell。"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


# 功能：运行 pySCENIC 官方命令，错误时立即停止并保留命令行信息。
def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


# 功能：读取主流程 h5ad、选取表达层并过滤极低表达基因。
def prepare_expression(args: argparse.Namespace) -> ad.AnnData:
    adata = ad.read_h5ad(args.input)
    if args.use_raw:
        if adata.raw is None:
            raise ValueError("--use-raw 已启用，但 h5ad 中没有 adata.raw。")
        expression = adata.raw.to_adata()
        expression.obs = adata.obs.copy()
    else:
        expression = adata.copy()
        if args.layer:
            expression.X = expression.layers[args.layer].copy()
    detected = np.asarray((expression.X > 0).sum(axis=0)).ravel()
    expression = expression[:, detected >= args.min_cells].copy()
    expression.obs_names_make_unique()
    expression.var_names_make_unique()
    # AnnData 默认写出 obs_names/var_names；pySCENIC CLI 默认读取 CellID/Gene。
    expression.obs_names.name = "CellID"
    expression.var_names.name = "Gene"
    return expression


# 功能：依次执行 GRNBoost2、cisTarget motif pruning 和 AUCell 三个官方步骤。
def run_pyscenic(args: argparse.Namespace, loom_path: Path) -> tuple[Path, Path, Path]:
    adjacency = args.output / "01_grnboost2_adjacencies.tsv"
    regulons = args.output / "02_ctx_regulons.csv"
    auc_loom = args.output / "03_aucell.loom"
    run_command(
        [
            "pyscenic", "grn", str(loom_path), str(args.tf_list),
            "--method", "grnboost2", "--output", str(adjacency),
            "--num_workers", str(args.workers), "--seed", str(args.seed),
        ]
    )
    # ranking 数据库与 motif annotations 必须属于同一物种和相容版本。
    run_command(
        [
            "pyscenic", "ctx", str(adjacency),
            *[str(path) for path in args.databases],
            "--annotations_fname", str(args.motif_annotations),
            "--expression_mtx_fname", str(loom_path),
            "--output", str(regulons), "--num_workers", str(args.workers),
        ]
    )
    run_command(
        [
            "pyscenic", "aucell", str(loom_path), str(regulons),
            "--output", str(auc_loom), "--num_workers", str(args.workers),
            "--seed", str(args.seed),
        ]
    )
    return adjacency, regulons, auc_loom


# 功能：从 AUCell loom 的结构化列属性提取 cells×regulons 活性矩阵。
def read_regulon_auc(auc_loom: Path) -> pd.DataFrame:
    import loompy

    with loompy.connect(str(auc_loom), mode="r") as dataset:
        auc = pd.DataFrame(
            dataset.ca.RegulonsAUC,
            index=pd.Index(dataset.ca.CellID.astype(str), name="cell"),
        )
    auc.columns = [str(column) for column in auc.columns]
    return auc


# 功能：按细胞类型计算 regulon specificity score。
def calculate_rss(auc: pd.DataFrame, obs: pd.DataFrame, celltype_col: str) -> pd.DataFrame:
    from pyscenic.rss import regulon_specificity_scores

    rss = regulon_specificity_scores(auc, obs.loc[auc.index, celltype_col].astype(str))
    rss.index.name = "cell_type"
    return rss


# 功能：同时保存 600 dpi PNG 和矢量 PDF。
def save_figure(figure: plt.Figure, stem: Path, width: float, height: float) -> None:
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(figure)


# 功能：在主流程固定 UMAP 坐标上展示方差最高的 regulon 活性。
def plot_top_regulons(adata: ad.AnnData, auc: pd.DataFrame, output: Path, top_n: int) -> None:
    if "X_umap" not in adata.obsm:
        return
    aligned = auc.reindex(adata.obs_names)
    top = aligned.var(axis=0).sort_values(ascending=False).head(top_n).index
    ncol = 3
    nrow = int(np.ceil(len(top) / ncol))
    figure, axes = plt.subplots(nrow, ncol, squeeze=False)
    coordinates = adata.obsm["X_umap"]
    cmap = LinearSegmentedColormap.from_list(
        "npg_sequential", ["#FFF5F0", "#F39B7F", "#E64B35", "#7A0019"]
    )
    for axis, regulon in zip(axes.ravel(), top, strict=False):
        points = axis.scatter(
            coordinates[:, 0], coordinates[:, 1], c=aligned[regulon],
            s=2, cmap=cmap, rasterized=False,
        )
        axis.set_title(regulon)
        axis.set_aspect("equal")
        axis.set_axis_off()
        figure.colorbar(points, ax=axis, fraction=0.04)
    for axis in axes.ravel()[len(top):]:
        axis.set_axis_off()
    save_figure(figure, output / "top_variable_regulon_umap", 12, 4 * nrow)


# 功能：定义 h5ad、TF 列表、cisTarget 数据库和运行参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run pySCENIC from an annotated h5ad.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tf-list", type=Path, required=True)
    parser.add_argument("--databases", type=Path, nargs="+", required=True)
    parser.add_argument("--motif-annotations", type=Path, required=True)
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use-raw", action="store_true")
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1001)
    parser.add_argument("--top-regulons", type=int, default=12)
    return parser


# 功能：串联 pySCENIC、RSS、h5ad 写回、绘图和中文结果说明。
def main() -> None:
    args = build_parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    expression = prepare_expression(args)
    loom_path = args.output / "00_expression.loom"
    expression.write_loom(loom_path)
    _, _, auc_loom = run_pyscenic(args, loom_path)
    auc = read_regulon_auc(auc_loom)
    auc.to_csv(args.output / "regulon_auc.csv", index=True)

    original = ad.read_h5ad(args.input)
    if args.celltype_col in original.obs:
        calculate_rss(auc, original.obs, args.celltype_col).to_csv(
            args.output / "regulon_specificity_scores.csv", index=True
        )
    original.obsm["X_regulon_auc"] = auc.reindex(original.obs_names).to_numpy()
    original.uns["regulon_names"] = auc.columns.to_list()
    original.write_h5ad(args.output / "pyscenic_annotated.h5ad")
    plot_top_regulons(original, auc, args.output, args.top_regulons)

    (args.output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# pySCENIC 结果说明", "",
                f"- 输入 h5ad：{args.input}",
                f"- 基因过滤：至少在 {args.min_cells} 个细胞中表达。",
                f"- TF 列表：{args.tf_list}",
                f"- cisTarget 数据库：{', '.join(map(str, args.databases))}",
                f"- motif annotations：{args.motif_annotations}",
                "- 01_grnboost2_adjacencies.tsv：TF-target 候选共表达边。",
                "- 02_ctx_regulons.csv：经 motif 富集和剪枝后的 regulon。",
                "- regulon_auc.csv：每个细胞的 AUCell regulon 活性。",
                "- regulon_specificity_scores.csv：按细胞类型计算的 RSS。",
                "- SCENIC 推断的是调控关联；motif 支持不等同于实验验证的 TF 直接结合。",
            ]
        ) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
