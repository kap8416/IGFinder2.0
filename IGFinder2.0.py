#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
IGFinder_final_debug.py

Safe version that analyzes only chromosomes 1 to 22 to avoid Ensembl API errors.
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
    parser = argparse.ArgumentParser(description='IGFinder debug: intronless genes with visible messages.')
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
    logging.info(f'🧬 🧬 Retrieving genes for {species} desde Ensembl...')
    chrom_data, _ = fetch_genes.chromosomes_info(species=species)

    # Mantener solo cromosomas 1 al 22
    valid_chroms = [str(i) for i in range(1, 23)]
    chrom_data = {k: v for k, v in chrom_data.items() if k in valid_chroms}

    all_genes = []
    for chrom in chrom_data:
        logging.info(f"→ → Processing chromosome {chrom}...")
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
    plt.title("Gene size by type")
    plt.savefig(f"{output_prefix}_boxplot.png"); plt.close()

    plt.figure(figsize=(8,5))
    sns.violinplot(x="type", y="length", data=df)
    plt.title("Gene length distribution")
    plt.savefig(f"{output_prefix}_violin.png"); plt.close()

    plt.figure(figsize=(8,5))
    sns.kdeplot(data=df[df['type'] == 'intronless']['length'], label="Intronless", fill=True)
    sns.kdeplot(data=df[df['type'] == 'multi-exonic']['length'], label="Multi-exonic", fill=True)
    plt.title("Gene size density curves")
    plt.legend(); plt.savefig(f"{output_prefix}_density.png"); plt.close()

    plt.figure(figsize=(10,5))
    sns.countplot(data=df, x='chr', hue='type')
    plt.title("Gene distribution per chromosome")
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
