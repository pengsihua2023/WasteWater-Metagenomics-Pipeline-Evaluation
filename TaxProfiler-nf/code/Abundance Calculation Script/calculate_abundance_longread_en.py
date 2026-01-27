#!/usr/bin/env python3
"""
Long-read data abundance calculation script (based on Kraken2 report)

For long-read data (Nanopore/PacBio), Kraken2 classification is already accurate enough,
and Bracken statistical correction is not needed. Extract abundance directly from Kraken2 report.

Usage:
    python calculate_abundance_longread_en.py \
        --kraken results/kraken2/sample.report \
        --output sample_abundance.tsv
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_kraken_report(kraken_file, taxonomy_level='S'):
    """
    Parse Kraken2 report file, extract abundance at specified taxonomic level
    
    Kraken2 report format:
    percentage  clade_reads  taxon_reads  rank  taxid  name
    """
    print(f"📖 Reading Kraken2 report: {kraken_file}")
    
    try:
        # Read Kraken2 report (no header, tab-separated)
        df = pd.read_csv(
            kraken_file, 
            sep='\t', 
            header=None,
            names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name']
        )
        
        # Clean species names (remove leading spaces and rank markers)
        df['name'] = df['name'].str.strip()
        
        # Filter to specified taxonomic level
        # S = Species, G = Genus, F = Family, O = Order, C = Class, P = Phylum, D = Domain
        df_filtered = df[df['rank'] == taxonomy_level].copy()
        
        print(f"✅ Found {len(df_filtered)} {get_level_name(taxonomy_level)}-level taxonomic units")
        
        return df_filtered
        
    except Exception as e:
        print(f"❌ Error: Unable to read Kraken2 report - {e}")
        sys.exit(1)


def get_level_name(level_code):
    """Get taxonomic level name"""
    level_map = {
        'D': 'domain',
        'P': 'phylum',
        'C': 'class',
        'O': 'order',
        'F': 'family',
        'G': 'genus',
        'S': 'species'
    }
    return level_map.get(level_code, level_code)


def get_genome_lengths():
    """
    Get genome lengths of common viruses (bp)
    Consistent with short-read script
    """
    genome_lengths = {
        # Coronaviruses
        'Severe acute respiratory syndrome coronavirus 2': 29903,
        'SARS-CoV-2': 29903,
        'Human coronavirus 229E': 27317,
        'Human coronavirus OC43': 30738,
        'Middle East respiratory syndrome-related coronavirus': 30119,
        
        # Influenza viruses
        'Influenza A virus': 13588,
        'Influenza B virus': 14548,
        
        # Herpesviruses
        'Human alphaherpesvirus 1': 152261,
        'Human alphaherpesvirus 2': 154746,
        'Human betaherpesvirus 5': 235646,
        'Epstein-Barr virus': 171823,
        
        # Adenoviruses
        'Human adenovirus C': 35935,
        
        # HIV
        'Human immunodeficiency virus 1': 9181,
        'Human immunodeficiency virus 2': 10359,
        
        # Hepatitis viruses
        'Hepatitis B virus': 3182,
        'Hepatitis C virus': 9646,
        
        # Other common viruses
        'Human papillomavirus': 7904,
        'Norwalk virus': 7654,
        'Rotavirus A': 18555,
        'Respiratory syncytial virus': 15222,
    }
    
    return genome_lengths


def calculate_abundance(kraken_df, genome_lengths):
    """
    Calculate RPM and RPKM
    
    For long-read data:
    - Use clade_reads from Kraken2 report (includes reads from sub-taxonomic units)
    - percentage is already relative abundance (0-100)
    """
    print("\n📊 Calculating abundance metrics...")
    
    # Calculate total reads (use root node reads)
    # Note: Kraken2 report percentages are relative to total reads
    # We need to back-calculate total reads from percentages
    if len(kraken_df) > 0:
        # Calculate total reads using clade_reads and percentage
        max_clade = kraken_df['clade_reads'].max()
        max_pct = kraken_df[kraken_df['clade_reads'] == max_clade]['percentage'].values[0]
        if max_pct > 0:
            total_reads = int(max_clade * 100 / max_pct)
        else:
            total_reads = kraken_df['clade_reads'].sum()
    else:
        print("⚠️  Warning: No classification data found")
        return pd.DataFrame()
    
    total_reads_million = total_reads / 1_000_000
    
    print(f"   Total reads (estimated): {total_reads:,}")
    print(f"   Total reads (million): {total_reads_million:.2f} M")
    
    # Create results dataframe
    results = []
    
    for idx, row in kraken_df.iterrows():
        species_name = row['name']
        taxonomy_id = row['taxid']
        assigned_reads = row['clade_reads']  # Use clade_reads
        fraction = row['percentage'] / 100  # Convert to 0-1
        
        # Calculate RPM
        rpm = (assigned_reads / total_reads) * 1_000_000
        
        # Try to get genome length
        genome_length = None
        rpkm = None
        
        # Try exact match
        if species_name in genome_lengths:
            genome_length = genome_lengths[species_name]
        else:
            # Try partial match
            for key in genome_lengths:
                if key in species_name or species_name in key:
                    genome_length = genome_lengths[key]
                    break
        
        # If genome length found, calculate RPKM
        if genome_length:
            genome_length_kb = genome_length / 1000
            rpkm = assigned_reads / (genome_length_kb * total_reads_million)
        
        results.append({
            'Species': species_name,
            'Taxonomy_ID': taxonomy_id,
            'Assigned_Reads': int(assigned_reads),
            'Fraction': fraction,
            'RPM': rpm,
            'Genome_Length_bp': genome_length if genome_length else 'NA',
            'RPKM': rpkm if rpkm else 'NA'
        })
    
    result_df = pd.DataFrame(results)
    
    # Sort by RPM in descending order
    result_df = result_df.sort_values('RPM', ascending=False)
    
    print(f"✅ Calculation complete, total {len(result_df)} species")
    
    return result_df


def main():
    parser = argparse.ArgumentParser(
        description='Long-read data abundance calculation (based on Kraken2 report)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Single sample
  python calculate_abundance_longread_en.py \\
    --kraken results/kraken2/sample1.report \\
    --output sample1_abundance.tsv
  
  # Specify taxonomic level
  python calculate_abundance_longread_en.py \\
    --kraken results/kraken2/sample1.report \\
    --level G \\
    --output sample1_abundance_genus.tsv

Taxonomic level codes:
  D = Domain
  P = Phylum
  C = Class
  O = Order
  F = Family
  G = Genus
  S = Species (default)
        """
    )
    
    parser.add_argument('--kraken', required=True,
                       help='Kraken2 report file path (*.report)')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--level', default='S',
                       help='Taxonomic level (default: S=species)')
    parser.add_argument('--genome-db', 
                       help='Custom genome length database (TSV: species<tab>length)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("🧬 Long-read Viral Abundance Calculator (RPM & RPKM)")
    print("   Based on Kraken2 report (no Bracken needed)")
    print("="*60)
    
    # Read Kraken2 report
    kraken_df = parse_kraken_report(args.kraken, args.level)
    
    if len(kraken_df) == 0:
        print("\n⚠️  Warning: No classification data found")
        # Create empty result file
        pd.DataFrame(columns=['Species', 'Taxonomy_ID', 'Assigned_Reads', 
                             'Fraction', 'RPM', 'Genome_Length_bp', 'RPKM']).to_csv(
            args.output, sep='\t', index=False
        )
        sys.exit(0)
    
    # Get genome lengths
    genome_lengths = get_genome_lengths()
    
    # If custom genome database provided, read and merge
    if args.genome_db:
        print(f"\n📖 Reading custom genome length database: {args.genome_db}")
        try:
            custom_db = pd.read_csv(args.genome_db, sep='\t', 
                                   names=['species', 'length'])
            custom_lengths = dict(zip(custom_db['species'], custom_db['length']))
            genome_lengths.update(custom_lengths)
            print(f"✅ Added {len(custom_lengths)} custom genome lengths")
        except Exception as e:
            print(f"⚠️  Warning: Unable to read custom database - {e}")
    
    # Calculate abundance
    result_df = calculate_abundance(kraken_df, genome_lengths)
    
    # Save results
    print(f"\n💾 Saving results to: {args.output}")
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # Print statistical summary
    print("\n" + "="*60)
    print("📈 Statistical Summary")
    print("="*60)
    print(f"Detected {get_level_name(args.level)} count: {len(result_df)}")
    print(f"{get_level_name(args.level).capitalize()} with genome length info: {(result_df['RPKM'] != 'NA').sum()}")
    print(f"\nTop 10 most abundant {get_level_name(args.level)} (by RPM):")
    print("-"*60)
    
    top10 = result_df.head(10)
    for idx, row in top10.iterrows():
        print(f"{row['Species'][:50]:50s} RPM: {row['RPM']:10.2f}")
    
    print("\n✅ Complete!")


if __name__ == '__main__':
    main()
