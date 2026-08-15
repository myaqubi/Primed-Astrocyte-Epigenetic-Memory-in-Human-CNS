from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


WORKBOOK = Path("Astrocyte_memory_genes_qPCR_data.xlsx")
OUTPUT_DIR = Path("qPCR_memory_heatmap")
CONDITIONS = ["Ctrl", "Stim", "Prime", "Dhit"]

# Source donor/sample IDs are the columns in the qPCR workbook.
SOURCE_DONORS = ["HFA598", "HFA596", "HFA593", "HFA579"]

# Display names match the final heatmap donor renaming.
DISPLAY_DONORS = {
    "HFA593": "Donor1",
    "HFA579": "Donor2",
    "HFA596": "Donor3",
    "HFA598": "Donor4",
}
DISPLAY_ORDER = ["HFA593", "HFA579", "HFA596", "HFA598"]

DONOR_STYLES = {
    "Donor1": {"color": "#000000", "marker": "D"},
    "Donor2": {"color": "#E69F00", "marker": "o"},
    "Donor3": {"color": "#56B4E9", "marker": "s"},
    "Donor4": {"color": "#CC79A7", "marker": "^"},
}


def safe_name(text: str) -> str:
    return "".join(char if char.isalnum() or char in "-_" else "_" for char in text)


def read_gene_blocks() -> list[dict]:
    raw = pd.read_excel(WORKBOOK, header=None)
    blocks = []

    for row_idx in range(len(raw)):
        gene = raw.iloc[row_idx, 1] if raw.shape[1] > 1 else None
        if pd.isna(gene):
            continue

        data_row_start = row_idx + 2
        data_row_end = row_idx + 6
        if data_row_end > len(raw):
            continue

        donor_row = raw.iloc[row_idx + 1, 1:5].tolist()
        if any(pd.isna(value) for value in donor_row):
            continue

        data = raw.iloc[data_row_start:data_row_end, 0:5].copy()
        data.columns = ["Condition", *SOURCE_DONORS]
        data = data[data["Condition"].isin(CONDITIONS)]
        if set(data["Condition"]) != set(CONDITIONS):
            continue

        data["Condition"] = pd.Categorical(data["Condition"], CONDITIONS, ordered=True)
        data = data.sort_values("Condition")
        blocks.append({"gene": str(gene).strip(), "data": data})

    return blocks


def plot_gene(block: dict) -> None:
    gene = block["gene"]
    data = block["data"]
    x_values = list(range(len(CONDITIONS)))

    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

    for source_donor in DISPLAY_ORDER:
        display_donor = DISPLAY_DONORS[source_donor]
        style = DONOR_STYLES[display_donor]
        ax.plot(
            x_values,
            data[source_donor].astype(float),
            color=style["color"],
            marker=style["marker"],
            linewidth=2.1,
            markersize=6.5,
        )

    ax.set_title(gene, color="black", fontsize=24, fontweight="bold", pad=10)
    ax.set_ylabel("Fold change", color="black", fontsize=21, fontweight="bold")
    ax.set_xlabel("")
    ax.set_xticks(x_values)
    ax.set_xticklabels(CONDITIONS, rotation=35, ha="right")
    ax.tick_params(axis="both", colors="black", labelsize=16, width=2.0, length=6)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("bottom", "left"):
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(2.0)

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    fig.subplots_adjust(left=0.28, right=0.96, bottom=0.25, top=0.82)
    fig.savefig(
        OUTPUT_DIR / f"{safe_name(gene)}_lineplot_renamed_donors.png",
        facecolor="white",
    )
    plt.close(fig)


def save_legend() -> None:
    handles = []
    labels = []
    for source_donor in DISPLAY_ORDER:
        display_donor = DISPLAY_DONORS[source_donor]
        style = DONOR_STYLES[display_donor]
        handles.append(
            Line2D(
                [0],
                [0],
                color=style["color"],
                marker=style["marker"],
                linewidth=2.8,
                markersize=9,
            )
        )
        labels.append(display_donor)

    fig, ax = plt.subplots(figsize=(2.8, 1.8), dpi=300)
    ax.axis("off")
    legend = ax.legend(
        handles,
        labels,
        loc="center",
        frameon=False,
        fontsize=16,
        handlelength=2.0,
        handletextpad=0.8,
    )
    for text in legend.get_texts():
        text.set_fontweight("bold")
        text.set_color("black")

    fig.patch.set_facecolor("white")
    fig.savefig(
        OUTPUT_DIR / "lineplot_renamed_donors_legend.png",
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    for block in read_gene_blocks():
        plot_gene(block)
        print(f"Saved {block['gene']}")
    save_legend()
    print("Saved legend")


if __name__ == "__main__":
    main()
