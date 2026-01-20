#!/bin/bash

# Quick start script for testing structural variant-feature relationships
# This script runs the knowledge graph creation with just the new SV-feature adapters

set -e  # Exit on error

echo "=================================================="
echo "  Structural Variant-Feature Relationships Test  "
echo "=================================================="
echo ""

# Check if output directory argument is provided
OUTPUT_DIR=${1:-"./output_sv_features_test"}

echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if dependencies are installed
if ! command -v python &> /dev/null; then
    echo "Error: Python is not installed"
    exit 1
fi

# Check if the required files exist
SAMPLE_CONFIG="config/hsa/hsa_adapters_config_sample.yaml"
if [ ! -f "$SAMPLE_CONFIG" ]; then
    echo "Error: Sample config not found at $SAMPLE_CONFIG"
    exit 1
fi

echo "Step 1: Validating configuration files..."
python -c "import yaml; yaml.safe_load(open('$SAMPLE_CONFIG'))" || {
    echo "Error: Invalid YAML in sample config"
    exit 1
}
echo "✓ Configuration files are valid"
echo ""

echo "Step 2: Checking sample data files..."
REQUIRED_FILES=(
    "samples/hsa/dbvar_sample.vcf.gz"
    "samples/hsa/dgv_GRCh38_hg38_variants.txt.gz"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Warning: Sample file not found: $file"
        echo "  The adapter will skip this source if file is missing"
    else
        echo "✓ Found: $file"
    fi
done
echo ""

echo "Step 3: Running coordinate logic tests..."
if [ -f "test_coordinate_logic.py" ]; then
    python test_coordinate_logic.py || {
        echo "Error: Coordinate logic tests failed"
        exit 1
    }
    echo ""
else
    echo "Warning: test_coordinate_logic.py not found, skipping tests"
    echo ""
fi

echo "Step 4: Creating output directory..."
mkdir -p "$OUTPUT_DIR"
echo "✓ Output directory created: $OUTPUT_DIR"
echo ""

echo "Step 5: Running knowledge graph creation..."
echo "  This will process only the SV-feature relationship adapters"
echo "  Estimated time: 2-5 minutes"
echo ""

python create_knowledge_graph.py \
    --species hsa \
    --dataset sample \
    --output-dir "$OUTPUT_DIR" \
    --writer-type metta \
    --include-adapters \
        feature_located_in_structural_variant \
        feature_overlaps_structural_variant \
    || {
        echo ""
        echo "Error: Knowledge graph creation failed"
        echo "Check the logs above for details"
        exit 1
    }

echo ""
echo "=================================================="
echo "                  SUCCESS!                        "
echo "=================================================="
echo ""
echo "Output files are in: $OUTPUT_DIR"
echo ""
echo "To view the results:"
echo "  1. Graph info:  cat $OUTPUT_DIR/graph_info.json | jq ."
echo "  2. List files:  ls -lh $OUTPUT_DIR/"
echo ""
echo "To run with all adapters (complete KG):"
echo "  python create_knowledge_graph.py \\"
echo "    --species hsa \\"
echo "    --dataset sample \\"
echo "    --output-dir ./output_complete \\"
echo "    --writer-type metta"
echo ""
echo "See docs/STRUCTURAL_VARIANT_FEATURES_GUIDE.md for more options"
echo ""
