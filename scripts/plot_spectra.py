import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

COLOR_MUTATION = {
    "C>A": "#E64B35", "C>G": "#4DBBD5", "C>T": "#00A087",
    "T>A": "#3C5488", "T>C": "#F39B7F", "T>G": "#8491B4"
}

COLOR_TRIPLET = {
    "C": "#7E6148", "T": "#B09C85"
}


def plot_mutations(series, output_path, title):
    df = series.reset_index()
    df.columns = ["index", "count"]
    df[['First_Base', 'Mutation', 'Third_Base']] = df['index'].str.extract(r'(\w)\[(\w>\w)\](\w)')
    df = df.sort_values(by=['Mutation', 'First_Base', 'Third_Base'])
    df.index = df["index"]
    sorted_data = df["count"]

    colors = [COLOR_MUTATION[m.split("[")[1][:-2]] for m in sorted_data.index]

    plt.figure(figsize=(12, 5), dpi=300)
    plt.bar(sorted_data.index, sorted_data.values, color=colors)
    plt.xticks(rotation=90, fontsize=6)
    plt.ylabel("Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_triplets(series, output_path, title):
    df = series.copy()
    df = df.sort_index()
    colors = [COLOR_TRIPLET.get(t[1], "gray") for t in df.index]

    plt.figure(figsize=(12, 5), dpi=300)
    plt.bar(df.index, df.values, color=colors)
    plt.xticks(rotation=90, fontsize=6)
    plt.ylabel("Triplet Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_mutation_spectra_overlay(data1, data2, labels, file_name=None):
    """
    Overlays mutation spectra for two datasets with color-coded x-axis categories.
    
    Parameters:
        data1: pd.Series - Mutation counts for species 1.
        data2: pd.Series - Mutation counts for species 2.
        labels: list - [label1, label2] names for species being compared.
        file_name: str - Optional. Path to save the figure.
    """
    def prepare_data(data):
        df = data.reset_index()
        df.columns = ["index", "count"]
        df[['First_Base', 'Mutation', 'Third_Base']] = df['index'].str.extract(r'(\w)\[(\w>\w)\](\w)')
        df = df.sort_values(by=['Mutation', 'First_Base', 'Third_Base'])
        df.index = df["index"]
        return df["count"]

    sorted_data1 = prepare_data(data1)
    sorted_data2 = prepare_data(data2)

    categories = sorted_data1.index
    color_dict = {"C>A": "red", "C>G": "green", "C>T": "blue", "T>A": "orange", "T>C": "purple", "T>G": "brown"}
    tick_colors = [color_dict[m.split("[")[1][:-2]] for m in categories]

    fig, ax = plt.subplots(figsize=(14, 6), dpi=300)
    x = range(len(categories))

    ax.bar(x, sorted_data2.values, color="red", width=0.6, label=labels[1], align='center', alpha=0.5)
    ax.bar(x, sorted_data1.values, color="yellow", width=0.6, label=labels[0], align='center', alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=6, rotation=90)
    for tick, color in zip(ax.get_xticklabels(), tick_colors):
        tick.set_color(color)

    plt.xlabel("Mutation category")
    plt.ylabel("Mutation count")
    plt.title(f"Mutation Spectra Comparison: {labels[0]} vs {labels[1]}")
    plt.legend()
    plt.tight_layout()
    if file_name:
        plt.savefig(file_name)
    else:
        plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot mutation and triplet spectra.")
    parser.add_argument("--input-dir", required=True, help="Directory containing the .tsv summary tables")
    parser.add_argument("--output-dir", help="Directory for the output plots")
    args = parser.parse_args()

    input_dir = args.input_dir
    if args.output_dir:
        plots_dir = args.output_dir
    else:
        plots_dir = os.path.join(os.path.dirname(input_dir), "Plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Load tables
    norm = pd.read_csv(os.path.join(input_dir, "normalized_scaled.tsv"), sep='\t', index_col=0)
    raw = pd.read_csv(os.path.join(input_dir, "collapsed_mutations.tsv"), sep='\t', index_col=0)
    scaled = pd.read_csv(os.path.join(input_dir, "scaled_raw.tsv"), sep='\t', index_col=0)
    trip = pd.read_csv(os.path.join(input_dir, "triplets.tsv"), sep='\t', index_col=0)

    for col in norm.columns:
        plot_mutations(norm[col], os.path.join(plots_dir, f"{col}_normalized.png"), f"Normalized Mutation Spectrum: {col}")

    for col in raw.columns:
        plot_mutations(raw[col], os.path.join(plots_dir, f"{col}_raw.png"), f"Raw Mutation Spectrum: {col}")

    for col in trip.columns:
        plot_triplets(trip[col], os.path.join(plots_dir, f"{col}_triplets.png"), f"Triplet Spectrum: {col}")

    # Optional: overlay comparison if there are exactly 2 species
    if len(norm.columns) == 2:
        species1, species2 = norm.columns
        plot_mutation_spectra_overlay(
            norm[species1], norm[species2],
            labels=[species1, species2],
            file_name=os.path.join(plots_dir, f"{species1}_vs_{species2}_normalized_overlay.png")
        )

        plot_mutation_spectra_overlay(
            scaled[species1], scaled[species2],
            labels=[species1, species2],
            file_name=os.path.join(plots_dir, f"{species1}_vs_{species2}_overlay.png")
        )



if __name__ == "__main__":
    main()