#!/bin/bash
# Simple abundance calculation script (using awk, no Python required)

INPUT=$1
OUTPUT=$2

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: ./calculate_abundance.sh <input CSV> <output file>"
    echo "Example: ./calculate_abundance.sh virus_classification.csv virus_abundance.csv"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: Input file does not exist: $INPUT"
    exit 1
fi

# Extract classification name column (assume column 3 is Classification_Name)
# Count occurrences of each classification and calculate abundance
awk -F',' '
BEGIN {
    # Skip CSV header
    OFS=","
    print "Classification,Read_Count,Relative_Abundance_%"
}
NR > 1 {
    # Assume Classification_Name is in column 3, adjust according to actual situation
    classification = $3
    # Remove possible quotes
    gsub(/^[" ]+|[" ]+$/, "", classification)
    if (classification != "" && classification != "Classification_Name") {
        count[classification]++
        total++
    }
}
END {
    # Calculate and output abundance
    for (class in count) {
        abundance = (count[class] / total) * 100
        printf "%s,%d,%.4f\n", class, count[class], abundance
    }
}' $INPUT | sort -t',' -k3 -rn > $OUTPUT

echo "Abundance calculation completed! Results saved in $OUTPUT"
echo "Total number of taxa: $(tail -n +2 $OUTPUT | wc -l)"

