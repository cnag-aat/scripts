#!/usr/bin/env python3
"""
get_closest_mitogenome.py

Given an NCBI taxon ID, walk up the taxonomy tree rank by rank
and find the closest available mitogenome(s) in RefSeq.
Stops at order. Returns up to 3 representative sequences.

Usage:
    python get_closest_mitogenome.py --taxid 9606 --email your@email.com
    python get_closest_mitogenome.py --taxid 9606 --email your@email.com --outdir ./output
"""

import argparse
import os
import sys
import time
from Bio import Entrez, SeqIO

# Ranks to walk through, from most specific to least (stopping at order)
RANKS_TO_SEARCH = ["species", "genus", "family", "order"]

MAX_RESULTS = 3


def fetch_lineage(taxid):
    """
    Fetch the full lineage for a taxid from NCBI Taxonomy.
    Returns a list of dicts with keys: TaxId, ScientificName, Rank
    ordered from most specific (species) to least (order/class/...).
    Also returns the name of the input taxon itself.
    """
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    if not records:
        print(f"ERROR: No taxonomy record found for taxid {taxid}", file=sys.stderr)
        sys.exit(1)

    record = records[0]
    input_name = record["ScientificName"]
    input_rank = record.get("Rank", "no rank")

    # LineageEx goes from root -> species, we want species -> order
    # Include the taxon itself at the front
    lineage = list(reversed(record["LineageEx"]))

    # Prepend the taxon itself so we search at species level first
    lineage.insert(0, {
        "TaxId": record["TaxId"],
        "ScientificName": input_name,
        "Rank": input_rank,
    })

    return lineage, input_name, input_rank


def search_mitogenomes_for_taxid(taxid, rank_name, sci_name):
    """
    Search nuccore for complete mitochondrial genomes under a given taxid.
    Returns a list of accession IDs (up to MAX_RESULTS).
    """
    query = (
        f"txid{taxid}[Organism:exp] "
        f"AND mitochondrion[Filter] "
        f"AND (complete genome[Title]" # OR complete sequence[Title])"
    )
    print(f"  Searching at {rank_name} ({sci_name}, txid{taxid})...")

    handle = Entrez.esearch(db="nuccore", term=query, retmax=100)
    result = Entrez.read(handle)
    handle.close()
    time.sleep(0.4)  # Be polite to NCBI

    count = int(result["Count"])
    ids = result["IdList"]

    if count == 0:
        print(f"    -> No mitogenomes found.")
        return []

    print(f"    -> {count} mitogenome(s) found. "
          f"{'Using all.' if count <= MAX_RESULTS else f'Selecting {MAX_RESULTS} representatives.'}")

    # If more than MAX_RESULTS, spread selection across the list
    # for rough diversity rather than just taking the first N
    if len(ids) > MAX_RESULTS:
        step = len(ids) // MAX_RESULTS
        ids = [ids[i * step] for i in range(MAX_RESULTS)]

    return ids


def fetch_and_save_sequences(ids, outdir, rank_name, sci_name):
    """
    Fetch FASTA sequences for a list of nuccore IDs and save each
    to a separate file in outdir.
    """
    os.makedirs(outdir, exist_ok=True)
    saved = []

    for uid in ids:
        print(f"  Fetching sequence for UID {uid}...")
        handle = Entrez.efetch(db="nuccore", id=uid, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        time.sleep(0.4)

        # Build a clean filename from the accession
        # FASTA id is like "NC_012920.1" — take the first word
        accession = record.id.split()[0].replace("/", "_")
        safe_rank = rank_name.replace(" ", "_")
        safe_name = sci_name.replace(" ", "_")
        filename = f"{accession}_{safe_rank}_{safe_name}.fasta"
        filepath = os.path.join(outdir, filename)

        SeqIO.write(record, filepath, "fasta")
        print(f"    -> Saved: {filepath}  ({len(record.seq)} bp)")
        saved.append(filepath)

    return saved


def main():
    parser = argparse.ArgumentParser(
        description="Fetch closest mitogenome(s) by taxonomy for a given NCBI taxon ID."
    )
    parser.add_argument("--taxid", required=True, help="NCBI taxon ID")
    parser.add_argument("--email", required=True, help="Email address for NCBI Entrez")
    parser.add_argument(
        "--outdir", default="./mito_closest", help="Output directory for FASTA files"
    )
    args = parser.parse_args()

    Entrez.email = args.email

    print(f"\nFetching taxonomy for taxid {args.taxid}...")
    lineage, input_name, input_rank = fetch_lineage(args.taxid)

    print(f"Input taxon: {input_name} (rank: {input_rank})")
    print(f"Lineage (most specific first): "
          + " > ".join(
              f"{n['ScientificName']} [{n['Rank']}]"
              for n in lineage
              if n["Rank"] in RANKS_TO_SEARCH
          ))
    print()

    # Walk up the tree, stopping at order
    for node in lineage:
        rank = node["Rank"]
        if rank not in RANKS_TO_SEARCH:
            continue

        ids = search_mitogenomes_for_taxid(node["TaxId"], rank, node["ScientificName"])

        if ids:
            print(f"\nFound mitogenomes at rank: {rank} ({node['ScientificName']})")
            saved = fetch_and_save_sequences(
                ids,
                outdir=args.outdir,
                rank_name=rank,
                sci_name=node["ScientificName"],
            )
            print(f"\nDone. {len(saved)} file(s) written to: {args.outdir}")
            return

    print(
        f"\nNo mitogenomes found for {input_name} up to and including order level.",
        file=sys.stderr,
    )
    sys.exit(1)


if __name__ == "__main__":
    main()