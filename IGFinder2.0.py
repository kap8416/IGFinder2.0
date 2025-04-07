#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
IGFinder_final_debug.py

Versión segura que analiza solo cromosomas 1 al 22 para evitar errores del API de Ensembl.
"""

import os
import sys
import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, chi2_contingency
import fetch_genes

def parse_arguments():
    parser = argparse.ArgumentParser(description='IGFinder debug: genes intronless con mensajes visibles.')
    parser.add_argument('--species', required=True)
    parser.add_argument('--utr_db', required=True)
    parser.add_argument('--output', default='genes_filtrados.tsv')
    parser.add_argument('--stats', action='store_true')
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--log', default='IGFinder_log.txt')
    return parser.parse_args()

def setup_logging(logfile):
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(sys.stdout)
        ]
    )

def is_intronless(transcript):
    return len(transcript.get('Exon', [])) == 1

def fetch_all_genes(species):
    logging.info(f'🧬 Recuperando genes para {species} desde Ensembl...')
    chrom_data, _ = fetch_genes.chromosomes_info(species=species)

    # Mantener solo cromosomas 1 al 22
    valid_chroms = [str(i) for i in range(1, 23)]
    chrom_data = {k: v for k, v in chrom_data.items() if k in valid_chroms}

    all_genes = []
    for chrom in chrom_data:
        logging.info(f"→ Procesando cromosoma {chrom}...")
        genes = fetch_genes.genes_in_chrom(species, chrom, chrom_data[chrom]['length'])
        logging.info(f"  - {chrom}: {len(genes)} genes")
        all_genes += genes
    logging.info(f'[✔] Total de genes recuperados: {len(all_genes)}')
    return all_genes

def classify_genes(gene_info, utr_introns):
    intronless, multiexonic = [], []
    for gene in gene_info:
        gene_id = gene['id']
        chr_ = gene['seq_region_name']
        start, end = gene['start'], gene['end']
        biotype = gene.get('biotype', 'NA')
        transcripts = gene.get('Transcript', [])

        if any(is_intronless(t) for t in transcripts):
            if gene_id not in utr_introns:
                intronless.append((gene_id, start, end, chr_, biotype, 'intronless'))
        else:
            multiexonic.append((gene_id, start, end, chr_, biotype, 'multi-exonic'))

    logging.info(f'[✔] Genes intronless válidos: {len(intronless)}')
    logging.info(f'[✔] Genes multiexónicos: {len(multiexonic)}')

    df = pd.DataFrame(intronless + multiexonic,
                      columns=['id', 'start', 'end', 'chr', 'biotype', 'type'])
    df['length'] = df['end'] - df['start']
    return df

def run_statistical_comparisons(df, output_prefix="IGFinder_stats"):
    df_intr = df[df["type"] == "intronless"]
    df_mult = df[df["type"] == "multi-exonic"]
    t_stat, p_val = ttest_ind(df_intr["length"], df_mult["length"], equal_var=False)
    cont_table = pd.crosstab(df["chr"], df["type"])
    chi2, chi2_p, _, _ = chi2_contingency(cont_table)
    with open(f"{output_prefix}.txt", "w") as f:
        f.write(f"T-test (longitud): t = {t_stat:.2f}, p = {p_val:.4e}\n")
        f.write(f"Chi-cuadrado (cromosomas): chi2 = {chi2:.2f}, p = {chi2_p:.4e}\n")
    logging.info(f"[✔] Estadísticas guardadas en {output_prefix}.txt")

def generate_visualizations(df, output_prefix="IGFinder_plots"):
    sns.set(style='whitegrid')

    plt.figure(figsize=(8,5))
    sns.boxplot(x="type", y="length", data=df)
    plt.title("Tamaño de genes por tipo")
    plt.savefig(f"{output_prefix}_boxplot.png"); plt.close()

    plt.figure(figsize=(8,5))
    sns.violinplot(x="type", y="length", data=df)
    plt.title("Distribución de longitud génica")
    plt.savefig(f"{output_prefix}_violin.png"); plt.close()

    plt.figure(figsize=(8,5))
    sns.kdeplot(data=df[df['type'] == 'intronless']['length'], label="Intronless", fill=True)
    sns.kdeplot(data=df[df['type'] == 'multi-exonic']['length'], label="Multi-exonic", fill=True)
    plt.title("Curvas de densidad del tamaño génico")
    plt.legend(); plt.savefig(f"{output_prefix}_density.png"); plt.close()

    plt.figure(figsize=(10,5))
    sns.countplot(data=df, x='chr', hue='type')
    plt.title("Distribución de genes por cromosoma")
    plt.xticks(rotation=90)
    plt.savefig(f"{output_prefix}_chromosomal_distribution.png"); plt.close()

    logging.info(f"[✔] Visualizaciones generadas como {output_prefix}_*.png")

def main():
    args = parse_arguments()
    setup_logging(args.log)

    if not os.path.isfile(args.utr_db):
        logging.error("Archivo UTR no encontrado.")
        sys.exit(1)

    try:
        utr_introns = pd.read_csv(args.utr_db, sep='\t', header=None)[2].dropna().unique()
        all_genes = fetch_all_genes(args.species)
        gene_info = fetch_genes.get_info(all_genes, fields=['id', 'start', 'end', 'seq_region_name', 'biotype', 'Transcript'])

        if not gene_info:
            logging.warning("No se recuperaron genes desde Ensembl.")
            return

        df = classify_genes(gene_info, utr_introns)
        if df.empty:
            logging.warning("No se generó ninguna entrada en el DataFrame final.")
            return

        df.to_csv(args.output, sep='\t', index=False)
        logging.info(f"[✔] Archivo final guardado: {args.output}")

        if args.stats:
            run_statistical_comparisons(df)
        if args.plots:
            generate_visualizations(df)

        logging.info("[🏁] Ejecución finalizada correctamente.")
    except Exception as e:
        logging.exception(f"[❌] Error durante la ejecución: {str(e)}")

if __name__ == '__main__':
    main()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, chi2_contingency

def save_results(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"[✔] Archivo guardado: {output_file}")

def generate_statistics(df, stats_file="IGFinder_stats.txt"):
    logging.info("[📉] Realizando pruebas estadísticas...")

    df_intronless = df[df["gene_type"] == "intronless"]
    df_multi = df[df["gene_type"] == "multi_exonic"]

    stats_text = ""

    # t-test for gene length
    t_stat, p_val = ttest_ind(df_intronless["length"], df_multi["length"], equal_var=False)
    stats_text += f"T-test (gene length):\nStatistic = {t_stat:.4f}, P-value = {p_val:.4e}\n"

    # chi-squared test for chromosomal distribution
    contingency = pd.crosstab(df["chromosome"], df["gene_type"])
    chi2, chi2_p, _, _ = chi2_contingency(contingency)
    stats_text += f"Chi-squared (chromosome distribution):\nStatistic = {chi2:.4f}, P-value = {chi2_p:.4e}\n"

    with open(stats_file, "w") as f:
        f.write(stats_text)

    logging.info(f"[✔] Estadísticas guardadas en {stats_file}")

def plot_gene_distributions(df, prefix="IGFinder_plots"):
    logging.info("[📊] Generando visualizaciones...")

    sns.set(style="whitegrid")

    # Boxplot
    plt.figure(figsize=(6, 4))
    sns.boxplot(x="gene_type", y="length", data=df)
    plt.title("Gene Length Distribution (Boxplot)")
    plt.savefig(f"{prefix}_boxplot.png")
    plt.close()

    # Violinplot
    plt.figure(figsize=(6, 4))
    sns.violinplot(x="gene_type", y="length", data=df)
    plt.title("Gene Length Distribution (Violin)")
    plt.savefig(f"{prefix}_violinplot.png")
    plt.close()

    # Density
    plt.figure(figsize=(6, 4))
    for label, group in df.groupby("gene_type"):
        sns.kdeplot(group["length"], label=label)
    plt.title("Gene Length Density")
    plt.legend()
    plt.savefig(f"{prefix}_density.png")
    plt.close()

    # Barplot by chromosome
    plt.figure(figsize=(10, 5))
    chrom_counts = df.groupby(["chromosome", "gene_type"]).size().unstack(fill_value=0)
    chrom_counts.plot(kind="bar", stacked=True)
    plt.title("Genes per Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{prefix}_barplot.png")
    plt.close()

    logging.info(f"[✔] Gráficas guardadas como {prefix}_*.png")

# Modificar main para guardar resultados y ejecutar análisis si se solicitan
def main():
    args = parse_arguments()
    setup_logging(args.log)

    all_genes = fetch_all_genes(args.species)

    filtered_df = classify_genes(all_genes, args.utr_db)
    logging.info(f"[✔] Total de genes recuperados: {len(filtered_df)}")

    if args.output:
        save_results(filtered_df, args.output)

    if args.stats:
        generate_statistics(filtered_df)

    if args.plots:
        plot_gene_distributions(filtered_df)

    logging.info("[🏁] Ejecución finalizada correctamente.")

if __name__ == "__main__":
    main()


def plot_advanced_gene_distributions(df, prefix="IGFinder_pub"):

    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.ticker as ticker

    sns.set(style="whitegrid", context="talk", font_scale=1.4)

    # Boxplot
    plt.figure(figsize=(8, 6))
    ax = sns.boxplot(x="gene_type", y="length", data=df, notch=True, linewidth=2.5, palette="Set2")
    ax.set_title("Gene Length Distribution", fontsize=18, weight='bold')
    ax.set_xlabel("Gene Type", fontsize=16)
    ax.set_ylabel("Gene Length (bp)", fontsize=16)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x/1000)}k"))
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{prefix}_boxplot_highres.png", dpi=300)
    plt.close()

    # Violinplot
    plt.figure(figsize=(8, 6))
    ax = sns.violinplot(x="gene_type", y="length", data=df, inner="quartile", linewidth=1.5, palette="Pastel1")
    ax.set_title("Gene Length Violin Plot", fontsize=18, weight='bold')
    ax.set_xlabel("Gene Type", fontsize=16)
    ax.set_ylabel("Gene Length (bp)", fontsize=16)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x/1000)}k"))
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{prefix}_violin_highres.png", dpi=300)
    plt.close()

    # Density
    plt.figure(figsize=(8, 6))
    for label, group in df.groupby("gene_type"):
        sns.kdeplot(group["length"], fill=True, label=label, alpha=0.5, linewidth=2)
    plt.title("Gene Length Density Curve", fontsize=18, weight='bold')
    plt.xlabel("Gene Length (bp)", fontsize=16)
    plt.ylabel("Density", fontsize=16)
    plt.legend(title="Gene Type")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{prefix}_density_highres.png", dpi=300)
    plt.close()

    # Barplot by chromosome
    plt.figure(figsize=(12, 6))
    chrom_counts = df.groupby(["chromosome", "gene_type"]).size().unstack(fill_value=0)
    chrom_counts = chrom_counts.loc[sorted(chrom_counts.index, key=lambda x: int(x) if x.isdigit() else x)]
    chrom_counts.plot(kind="bar", stacked=True, colormap="tab20", edgecolor='black')
    plt.title("Gene Distribution per Chromosome", fontsize=18, weight='bold')
    plt.xlabel("Chromosome", fontsize=16)
    plt.ylabel("Gene Count", fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{prefix}_chromosome_barplot_highres.png", dpi=300)
    plt.close()

    logging.info("[✔] Visualizaciones de publicación generadas.")
