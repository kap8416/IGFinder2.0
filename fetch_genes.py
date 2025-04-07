#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fetch_genes.py

Versión robusta con consultas en ventanas para evitar errores 400 por tamaño de cromosoma.
"""

from itertools import islice
import requests
import time
import logging

ENSEMBL_REST_API = "https://rest.ensembl.org"
WINDOW_SIZE = 1000000  # 10 millones de bases

def head(iterable, n=10):
    iterator = iter(iterable)
    yield from islice(iterator, n)

def chunks_of(iterable, chunk_size=10):
    iterator = iter(iterable)
    while True:
        next_ = list(head(iterator, chunk_size))
        if next_:
            yield next_
        else:
            break

def ensembl_get(endpoint, params=None):
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "IGFinder/1.0 (katia.avi@gmail.com)"
    }
    url = f"{ENSEMBL_REST_API}{endpoint}"
    response = requests.get(url, headers=headers, params=params)
    if not response.ok:
        logging.error(f"[❌] Falló la llamada a Ensembl: {url}")
        response.raise_for_status()
    return response.json()

def overlapping_features(species, region, feature='gene'):
    endpoint = f"/overlap/region/{species}/{region}"
    return ensembl_get(endpoint, {"feature": feature})

def get_lookup_batch(ids, expand=True):
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "IGFinder/1.0 (katia.avi@gmail.com)"
    }
    url = f"{ENSEMBL_REST_API}/lookup/id"
    response = requests.post(url, headers=headers, json={"ids": ids, "params": {"expand": expand}})
    if not response.ok:
        logging.warning("Fallo en batch lookup.")
        return {}
    return response.json()

def get_info(genes, fields=None):
    if fields is None:
        fields = ['id', 'start', 'end', 'seq_region_name', 'biotype', 'Transcript']
    all_info = []
    for chunk in chunks_of(genes, chunk_size=50):
        ids = [g["id"] for g in chunk]
        try:
            data = get_lookup_batch(ids)
            all_info.extend(data.values())
        except Exception as e:
            logging.warning(f"Error al obtener batch info: {e}")
            time.sleep(2)
    return all_info

def genes_in_chrom(species, chrom, length, feature='gene'):
    all_genes = []
    for start in range(1, length, WINDOW_SIZE):
        end = min(start + WINDOW_SIZE - 1, length)
        region = f"{chrom}:{start}-{end}"
        try:
            genes = overlapping_features(species, region, feature)
            all_genes.extend(genes)
        except Exception as e:
            logging.warning(f"Fallo en ventana {region}: {e}")
            continue
    return all_genes

def chromosomes_info(species):
    logging.info(f"→ Conectando a Ensembl para info de ensamblado de {species}...")
    response = ensembl_get(f"/info/assembly/{species}")
    regions = response.get("top_level_region", [])
    chrom_lengths = {
        region["name"]: {
            "length": region["length"],
            "coord_system": region.get("coord_system", "chromosome")
        }
        for region in regions if region.get("coord_system") == "chromosome"
    }
    return chrom_lengths, response
