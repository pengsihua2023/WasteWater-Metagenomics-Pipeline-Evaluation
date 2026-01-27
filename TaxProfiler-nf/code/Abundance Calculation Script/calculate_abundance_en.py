#!/usr/bin/env python3
"""
Calculate viral abundance metrics: RPM and RPKM

RPM (Reads Per Million): Number of reads assigned to the species per million reads
RPKM (Reads Per Kilobase Million): Normalized abundance accounting for genome length

Usage:
    python calculate_abundance_en.py --bracken results/bracken/sample_bracken.tsv \
                                      --kraken results/kraken2/sample.report \
                                      --output results/abundance_RPM_RPKM.tsv
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_bracken_output(bracken_file):
    """
    Parse Bracken output file
    
    Bracken output format:
    name    taxonomy_id    taxonomy_lvl    kraken_assigned_reads    added_reads    new_est_reads    fraction_total_reads
    """
    print(f"📖 Reading Bracken file: {bracken_file}")
    
    try:
        df = pd.read_csv(bracken_file, sep='\t')
        print(f"✅ Successfully read {len(df)} taxonomic units")
        return df
    except Exception as e:
        print(f"❌ Error: Unable to read Bracken file - {e}")
        sys.exit(1)


def parse_kraken_report(kraken_file):
    """
    Parse Kraken2 report file for genome length information
    
    Kraken2 report format (tab-separated):
    percentage    num_covered    num_assigned    rank    taxid    name
    """
    print(f"📖 Reading Kraken2 report: {kraken_file}")
    
    try:
        # Kraken2 report has no header
        df = pd.read_csv(kraken_file, sep='\t', header=None,
                        names=['percentage', 'num_covered', 'num_assigned', 
                               'rank', 'taxid', 'name'])
        # Clean species names (remove leading spaces and rank markers)
        df['name'] = df['name'].str.strip()
        print(f"✅ Successfully read Kraken2 report with {len(df)} entries")
        return df
    except Exception as e:
        print(f"❌ Error: Unable to read Kraken2 report - {e}")
        sys.exit(1)


def get_genome_lengths():
    """
    Get genome lengths of common viruses (bp)
    
    Note: This is an example dictionary. Update based on your viral database
    or automatically extract from database taxonomy files
    """
    # Common viral genome lengths (unit: bp)
    genome_lengths = {
        # Coronaviruses
        'Severe acute respiratory syndrome coronavirus 2': 29903,
        'SARS-CoV-2': 29903,
        'Human coronavirus 229E': 27317,
        'Human coronavirus OC43': 30738,
        'Middle East respiratory syndrome-related coronavirus': 30119,
        
        # Influenza viruses
        'Influenza A virus': 13588,  # Average length
        'Influenza B virus': 14548,
        
        # Herpesviruses
        'Human alphaherpesvirus 1': 152261,  # HSV-1
        'Human alphaherpesvirus 2': 154746,  # HSV-2
        'Human betaherpesvirus 5': 235646,   # CMV
        'Epstein-Barr virus': 171823,        # EBV
        
        # Adenoviruses
        'Human adenovirus C': 35935,
        
        # HIV
        'Human immunodeficiency virus 1': 9181,
        'Human immunodeficiency virus 2': 10359,
        
        # Hepatitis viruses
        'Hepatitis B virus': 3182,
        'Hepatitis C virus': 9646,
        
        # Other common viruses
        'Human papillomavirus': 7904,  # Average
        'Norwalk virus': 7654,
        'Rotavirus A': 18555,
        'Respiratory syncytial virus': 15222,
    }
    
    return genome_lengths


def calculate_rpm_rpkm(bracken_df, kraken_df, genome_lengths):
    """
    Calculate RPM and RPKM values
    
    RPM = (assigned_reads / total_reads) × 1,000,000
    RPKM = (assigned_reads / (genome_length_kb × total_reads_million))
    """
    print("\n📊 Calculating abundance metrics...")
    
    # Calculate total reads
    total_reads = bracken_df['new_est_reads'].sum()
    total_reads_million = total_reads / 1_000_000
    
    print(f"   Total reads: {total_reads:,}")
    print(f"   Total reads (million): {total_reads_million:.2f} M")
    
    # Create results dataframe
    results = []
    
    for idx, row in bracken_df.iterrows():
        species_name = row['name']
        taxonomy_id = row['taxonomy_id']
        assigned_reads = row['new_est_reads']
        fraction = row['fraction_total_reads']
        
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
        description='Calculate viral abundance metrics (RPM and RPKM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Single sample
  python calculate_abundance_en.py \\
    --bracken results/bracken/sample1_bracken.tsv \\
    --kraken results/kraken2/sample1.report \\
    --output sample1_abundance.tsv
  
  # Batch processing
  for sample in sample1 sample2 sample3; do
    python calculate_abundance_en.py \\
      --bracken results/bracken/${sample}_bracken.tsv \\
      --kraken results/kraken2/${sample}.report \\
      --output ${sample}_abundance.tsv
  done
        """
    )
    
    parser.add_argument('--bracken', required=True,
                       help='Bracken output file path (*.tsv)')
    parser.add_argument('--kraken', required=True,
                       help='Kraken2 report file path (*.report)')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--genome-db', 
                       help='Custom genome length database (TSV: species<tab>length)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("🧬 Viral Abundance Calculator (RPM & RPKM)")
    print("="*60)
    
    # Read input files
    bracken_df = parse_bracken_output(args.bracken)
    kraken_df = parse_kraken_report(args.kraken)
    
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
    result_df = calculate_rpm_rpkm(bracken_df, kraken_df, genome_lengths)
    
    # Save results
    print(f"\n💾 Saving results to: {args.output}")
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # Print statistical summary
    print("\n" + "="*60)
    print("📈 Statistical Summary")
    print("="*60)
    print(f"Detected species count: {len(result_df)}")
    print(f"Species with genome length info: {(result_df['RPKM'] != 'NA').sum()}")
    print(f"\nTop 10 most abundant species (by RPM):")
    print("-"*60)
    
    top10 = result_df.head(10)
    for idx, row in top10.iterrows():
        print(f"{row['Species'][:50]:50s} RPM: {row['RPM']:10.2f}")
    
    print("\n✅ Complete!")


if __name__ == '__main__':
    main()
