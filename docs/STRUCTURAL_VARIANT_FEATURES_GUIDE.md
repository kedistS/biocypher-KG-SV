# Running the Structural Variant-Feature Relationship Adapter

This guide explains how to run the BioCypher knowledge graph creation with the new structural variant-feature relationship adapters.

## Prerequisites

1. **Install dependencies**:
   ```bash
   # Install UV package manager (if not installed)
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Install project dependencies
   cd /home/user/biocypher-KG-SV
   uv sync
   ```

2. **Verify sample data exists**:
   ```bash
   ls -la samples/hsa/dgv_GRCh38_hg38_variants.txt.gz
   ls -la samples/hsa/dbvar_sample.vcf.gz
   ```

## Running Options

### Option 1: Using Species Mode (Recommended)

Run with the HSA species configuration using sample data:

```bash
# Run with sample dataset
python create_knowledge_graph.py \
  --species hsa \
  --dataset sample \
  --output-dir ./output_sv_features \
  --writer-type metta \
  --include-adapters feature_located_in_structural_variant feature_overlaps_structural_variant
```

**Explanation**:
- `--species hsa`: Use human (Homo sapiens) configuration
- `--dataset sample`: Use sample data files from `samples/hsa/`
- `--output-dir`: Where to write the knowledge graph files
- `--writer-type`: Output format (metta, neo4j, prolog, parquet, kgx, networkx)
- `--include-adapters`: Only run the specified adapters (optional)

### Option 2: Run All Adapters with SV Relationships

To generate a complete knowledge graph including the new SV-feature relationships:

```bash
python create_knowledge_graph.py \
  --species hsa \
  --dataset sample \
  --output-dir ./output_complete_kg \
  --writer-type metta
```

This will run ALL configured adapters including:
- Gene, transcript, exon nodes
- Structural variant nodes (dbVar, DGV)
- Feature-SV relationships (located_in, overlaps_with)
- All other configured relationships

### Option 3: Manual Mode (Full Dataset)

For production use with full datasets:

```bash
python create_knowledge_graph.py \
  --output-dir ./output_production \
  --adapters-config config/hsa/hsa_adapters_config.yaml \
  --dbsnp-rsids aux_files/hsa/dbsnp_rsids.pkl \
  --dbsnp-pos aux_files/hsa/dbsnp_pos.pkl \
  --schema-config config/hsa/hsa_schema_config.yaml \
  --writer-type parquet \
  --include-adapters \
    dgv_variant \
    dbvar_variant \
    gencode_gene \
    gencode_transcripts \
    gencode_exon \
    epd_promoter \
    rna_central_non_coding_rna \
    feature_located_in_structural_variant \
    feature_overlaps_structural_variant
```

**Note**: This requires full data files at `/mnt/hdd_2/abdu/biocypher_data/`

## Testing the New Adapters

### Quick Test with Sample Data

1. **Verify adapters are configured**:
   ```bash
   grep -A 20 "feature_located_in_structural_variant" config/hsa/hsa_adapters_config_sample.yaml
   ```

2. **Run with just the new adapters**:
   ```bash
   python create_knowledge_graph.py \
     --species hsa \
     --dataset sample \
     --output-dir ./test_sv_relationships \
     --writer-type metta \
     --include-adapters \
       feature_located_in_structural_variant \
       feature_overlaps_structural_variant
   ```

3. **Check the output**:
   ```bash
   ls -lh test_sv_relationships/
   cat test_sv_relationships/graph_info.json
   ```

## Configuration Files

The new adapters are configured in two places:

### 1. Adapter Configuration
**File**: `config/hsa/hsa_adapters_config.yaml` (or `hsa_adapters_config_sample.yaml` for samples)

```yaml
feature_located_in_structural_variant:
  adapter:
    module: biocypher_metta.adapters.hsa.structural_variant_feature_relationship_adapter
    cls: StructuralVariantFeatureRelationshipAdapter
    args:
      sv_sources:
        dbvar: /path/to/dbvar.vcf.gz
        dgv: /path/to/dgv_variants.txt.gz
      feature_sources:
        gene:
          filepath: /path/to/gencode.annotation.gtf.gz
          format: gtf
        # ... more feature types
      label: located_in
      taxon_id: 9606
  outdir: structural_variant_relationships
  nodes: False
  edges: True
```

### 2. Schema Configuration
**File**: `config/hsa/hsa_schema_config.yaml`

Defines the edge types:
- `feature located in structural variant` (biolink:located_in)
- `feature overlaps structural variant` (biolink:overlaps)

## Expected Output

The adapters will generate edges with these properties:

```python
{
  "chr": "chr1",
  "feature_start": 1000,
  "feature_end": 2000,
  "sv_start": 500,
  "sv_end": 3000,
  "source": "BioCypher-KG-SV",
  "sv_source": "dbVar"  # or "DGV"
}
```

### Example Edge Format (MeTTa)

```scheme
; Feature located_in structural variant
(located_in
  (gene "ENSEMBL:ENSG00000223972")
  (structural_variant "SO:0001059_chr1_500_3000_GRCh38")
  (props (chr "chr1") (feature_start 1000) (feature_end 2000)))

; Feature overlaps_with structural variant
(overlaps_with
  (gene "ENSEMBL:ENSG00000227232")
  (structural_variant "SO:0001059_dgv1n82")
  (props (chr "chr1") (feature_start 14500) (feature_end 29500)))
```

## Troubleshooting

### 1. Module not found errors
```bash
# Install dependencies
uv sync

# Or with pip
pip install -r requirements.txt
```

### 2. Sample data not found
```bash
# Verify sample data exists
ls samples/hsa/dgv_GRCh38_hg38_variants.txt.gz
ls samples/hsa/dbvar_sample.vcf.gz

# If missing, check sample adapters config
cat config/hsa/hsa_adapters_config_sample.yaml | grep -A 10 "dgv_variant"
```

### 3. No edges generated
Check logs for:
- "Loaded structural variants from N chromosomes"
- "Processing [feature_type] features from [filepath]"
- "Generated N [relationship_type] edges"

### 4. Memory issues with large datasets
Use the Parquet writer for better memory efficiency:
```bash
python create_knowledge_graph.py \
  --species hsa \
  --dataset full \
  --output-dir ./output \
  --writer-type parquet \
  --buffer-size 50000
```

## Performance Notes

- **Structural variants** are loaded into memory at initialization (organized by chromosome)
- **Features** are processed in streaming fashion to minimize memory usage
- **Expected processing time** (sample data): 2-10 minutes depending on system
- **Expected processing time** (full data): 30 minutes to several hours

## Validation

To validate the implementation:

```bash
# Run coordinate logic tests
python test_coordinate_logic.py

# Expected output:
# ============================================================
# Running Coordinate Comparison Logic Tests
# ============================================================
# Testing coordinate comparison logic...
# ✓ All _is_fully_contained tests passed
# ✓ All _has_overlap_but_not_contained tests passed
#
# Testing edge cases...
# ✓ All edge case tests passed
#
# ============================================================
# Test Results: 2 passed, 0 failed
# ============================================================
```

## Next Steps

After successful generation:

1. **Review the output**:
   ```bash
   cat output_sv_features/graph_info.json | jq .
   ```

2. **Import into your target system** (Neo4j, Prolog, etc.)

3. **Query the relationships**:
   - Find all genes located in structural variants
   - Find all features that overlap with structural variants
   - Analyze variant-gene associations

## Support

- **Documentation**: See main [README.md](README.md)
- **Issues**: Report at https://github.com/kedistS/biocypher-KG-SV/issues
- **Tests**: Run `python test_coordinate_logic.py` to verify logic

## Architecture Summary

```
Input Sources:
├── Structural Variants
│   ├── dbVar (VCF format)
│   └── DGV (TSV format)
└── Genomic Features
    ├── GENCODE (genes, transcripts, exons - GTF)
    ├── EPD (promoters - BED)
    └── RNAcentral (non-coding RNAs - BED)

Processing:
├── Load all SVs into memory (organized by chromosome)
├── Stream through each feature type
├── Compare coordinates (1-based, inclusive)
└── Generate edges based on overlap type

Output:
├── located_in edges (full containment)
└── overlaps_with edges (partial overlap)
```

---

**Implementation complete and ready for use!** ✓
