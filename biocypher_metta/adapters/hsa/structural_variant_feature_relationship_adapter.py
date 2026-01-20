"""
Adapter for creating inferred relationships between structural variants and genomic features
based on genomic coordinate overlap analysis.

This adapter creates two types of relationships:
1. located_in: Feature is fully contained within a structural variant
2. overlaps_with: Feature and structural variant share genomic coordinates but feature is not fully contained
"""

import gzip
from collections import defaultdict
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_regulatory_region_id


class StructuralVariantFeatureRelationshipAdapter(Adapter):
    """
    Creates relationships between genomic features and structural variants based on coordinate overlap.

    This adapter reads structural variant data from dbVar and DGV sources, and genomic feature data
    from various sources (GENCODE, RNAcentral, etc.), then creates edges based on coordinate comparisons.
    """

    # Index mappings for different file formats
    GTF_INDEX = {'chr': 0, 'feature_type': 2, 'coord_start': 3, 'coord_end': 4, 'info': 8}
    BED_INDEX = {'chr': 0, 'coord_start': 1, 'coord_end': 2}
    DBVAR_INDEX = {'chr': 0, 'coord_start': 1, 'id': 2, 'type': 4, 'info': 7}
    DGV_INDEX = {'variant_accession': 0, 'chr': 1, 'coord_start': 2, 'coord_end': 3, 'type': 5}

    # Supported variant types from dbVar
    DBVAR_VARIANT_TYPES = {'<CNV>', '<DEL>', '<DUP>', '<INS>', '<INV>'}

    def __init__(
        self,
        sv_sources,
        feature_sources,
        label,
        write_properties=True,
        add_provenance=True,
        taxon_id=9606
    ):
        """
        Initialize the adapter.

        Args:
            sv_sources: Dict with 'dbvar' and/or 'dgv' keys mapping to file paths
            feature_sources: Dict with feature type keys ('gene', 'transcript', 'exon', 'promoter', 'non_coding_rna')
                           mapping to file paths and formats
            label: The edge label ('located_in' or 'overlaps_with')
            write_properties: Whether to include edge properties
            add_provenance: Whether to include provenance information
            taxon_id: Taxonomy ID (default: 9606 for human)
        """
        self.sv_sources = sv_sources or {}
        self.feature_sources = feature_sources or {}
        self.label = label
        self.taxon_id = taxon_id
        self.source = 'BioCypher-KG-SV'

        # Storage for structural variants organized by chromosome for efficient lookup
        self.structural_variants = defaultdict(list)  # chr -> list of (id, start, end, source)

        super(StructuralVariantFeatureRelationshipAdapter, self).__init__(write_properties, add_provenance)

        # Load all structural variants into memory
        self._load_structural_variants()

    def _load_structural_variants(self):
        """Load all structural variants from configured sources into memory."""

        # Load dbVar variants
        if 'dbvar' in self.sv_sources:
            self._load_dbvar_variants(self.sv_sources['dbvar'])

        # Load DGV variants
        if 'dgv' in self.sv_sources:
            self._load_dgv_variants(self.sv_sources['dgv'])

        print(f"Loaded structural variants from {len(self.structural_variants)} chromosomes")
        for chr_name, variants in self.structural_variants.items():
            print(f"  {chr_name}: {len(variants)} variants")

    def _load_dbvar_variants(self, filepath):
        """Load structural variants from dbVar VCF file."""
        try:
            with gzip.open(filepath, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    data = line.strip().split('\t')
                    if len(data) < 8:
                        continue

                    variant_type_key = data[self.DBVAR_INDEX['type']]
                    if variant_type_key not in self.DBVAR_VARIANT_TYPES:
                        continue

                    chr_raw = data[self.DBVAR_INDEX['chr']]
                    chr_name = f'chr{chr_raw}'
                    start = int(data[self.DBVAR_INDEX['coord_start']])

                    # Extract END position from INFO field
                    info = data[self.DBVAR_INDEX['info']].split(';')
                    end = start
                    for item in info:
                        if item.startswith('END='):
                            end = int(item.split('=')[1])
                            break

                    # Create variant ID matching dbVar adapter format
                    variant_id = f"SO:0001059_{data[self.DBVAR_INDEX['id']]}"

                    self.structural_variants[chr_name].append({
                        'id': variant_id,
                        'start': start,
                        'end': end,
                        'source': 'dbVar'
                    })

        except Exception as e:
            print(f"Error loading dbVar variants: {e}")

    def _load_dgv_variants(self, filepath):
        """Load structural variants from DGV file."""
        try:
            with gzip.open(filepath, 'rt') as f:
                next(f)  # Skip header

                for line in f:
                    data = line.strip().split('\t')
                    if len(data) < 8:
                        continue

                    chr_raw = data[self.DGV_INDEX['chr']]
                    chr_name = f'chr{chr_raw}'
                    # DGV uses 0-based coordinates, convert to 1-based
                    start = int(data[self.DGV_INDEX['coord_start']]) + 1
                    end = int(data[self.DGV_INDEX['coord_end']])

                    # Create variant ID matching DGV adapter format
                    region_id = build_regulatory_region_id(chr_name, start, end)
                    variant_id = f"SO:0001059_{region_id}"

                    self.structural_variants[chr_name].append({
                        'id': variant_id,
                        'start': start,
                        'end': end,
                        'source': 'DGV'
                    })

        except Exception as e:
            print(f"Error loading DGV variants: {e}")

    @staticmethod
    def _is_fully_contained(feat_chr, feat_start, feat_end, sv_chr, sv_start, sv_end):
        """
        Check if a feature is fully contained within a structural variant.

        Args:
            feat_chr: Feature chromosome
            feat_start: Feature start position (1-based, inclusive)
            feat_end: Feature end position (1-based, inclusive)
            sv_chr: Structural variant chromosome
            sv_start: Structural variant start position (1-based, inclusive)
            sv_end: Structural variant end position (1-based, inclusive)

        Returns:
            True if feature is fully contained within the structural variant
        """
        if feat_chr != sv_chr:
            return False
        return feat_start >= sv_start and feat_end <= sv_end

    @staticmethod
    def _has_overlap_but_not_contained(feat_chr, feat_start, feat_end, sv_chr, sv_start, sv_end):
        """
        Check if a feature overlaps with a structural variant but is not fully contained.

        This includes cases where:
        - The structural variant is fully contained within the feature
        - The overlap is partial (feature extends beyond SV boundaries on one or both sides)

        Args:
            feat_chr: Feature chromosome
            feat_start: Feature start position (1-based, inclusive)
            feat_end: Feature end position (1-based, inclusive)
            sv_chr: Structural variant chromosome
            sv_start: Structural variant start position (1-based, inclusive)
            sv_end: Structural variant end position (1-based, inclusive)

        Returns:
            True if feature overlaps with SV but is not fully contained
        """
        if feat_chr != sv_chr:
            return False

        # Check for any overlap: intervals overlap if they don't end before the other starts
        has_overlap = not (feat_end < sv_start or feat_start > sv_end)

        if not has_overlap:
            return False

        # Check if fully contained
        is_contained = feat_start >= sv_start and feat_end <= sv_end

        # Return true only if overlaps but not fully contained
        return not is_contained

    def _generate_feature_id(self, feature_type, info_dict, raw_id=None):
        """
        Generate feature ID matching the format used by feature adapters.

        Args:
            feature_type: Type of feature ('gene', 'transcript', 'exon', 'promoter', 'non_coding_rna')
            info_dict: Dictionary of attributes from GTF/GFF
            raw_id: Raw ID (if already extracted)

        Returns:
            Formatted feature ID
        """
        if feature_type in ['gene', 'transcript', 'exon']:
            # GENCODE features use ENSEMBL IDs
            if raw_id is None:
                if feature_type == 'gene':
                    raw_id = info_dict.get('gene_id', '').split('.')[0]
                elif feature_type == 'transcript':
                    raw_id = info_dict.get('transcript_id', '').split('.')[0]
                elif feature_type == 'exon':
                    raw_id = info_dict.get('exon_id', '').split('.')[0]

            if not raw_id:
                return None

            # Handle PAR_Y genes
            if raw_id.endswith('_PAR_Y'):
                return f"ENSEMBL:{raw_id}"
            else:
                return f"ENSEMBL:{raw_id.split('.')[0]}"

        elif feature_type == 'promoter':
            # Promoters use SO:0000167 prefix with coordinate-based ID
            return raw_id  # Already formatted by caller

        elif feature_type == 'non_coding_rna':
            # Non-coding RNAs use RNAcentral IDs
            return raw_id  # Already formatted by caller

        return None

    def _parse_gtf_info(self, info_string):
        """Parse GTF/GFF info field into a dictionary."""
        info_dict = {}
        for item in info_string.strip().split(';'):
            item = item.strip()
            if not item:
                continue
            parts = item.split(' ', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().strip('"')
                info_dict[key] = value
        return info_dict

    def _process_gtf_features(self, filepath, feature_type):
        """
        Process features from a GENCODE GTF file.

        Args:
            filepath: Path to GTF file
            feature_type: Type of feature to extract ('gene', 'transcript', 'exon')

        Yields:
            Tuples of (feature_id, chr, start, end)
        """
        try:
            if filepath.endswith('.gz'):
                file_handle = gzip.open(filepath, 'rt')
            else:
                file_handle = open(filepath, 'rt')

            for line in file_handle:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                # Check if this is the feature type we're looking for
                line_feature_type = fields[self.GTF_INDEX['feature_type']]
                if line_feature_type != feature_type:
                    continue

                chr_name = fields[self.GTF_INDEX['chr']]
                start = int(fields[self.GTF_INDEX['coord_start']])
                end = int(fields[self.GTF_INDEX['coord_end']])

                # Parse info field
                info = self._parse_gtf_info(fields[self.GTF_INDEX['info']])

                # Generate feature ID
                feature_id = self._generate_feature_id(feature_type, info)

                if feature_id:
                    yield (feature_id, chr_name, start, end)

            file_handle.close()

        except Exception as e:
            print(f"Error processing GTF file for {feature_type}: {e}")

    def _process_bed_features(self, filepath, feature_type, id_column=3):
        """
        Process features from a BED format file.

        Args:
            filepath: Path to BED file
            feature_type: Type of feature ('promoter', 'non_coding_rna')
            id_column: Column index for feature ID (0-based)

        Yields:
            Tuples of (feature_id, chr, start, end)
        """
        try:
            if filepath.endswith('.gz'):
                file_handle = gzip.open(filepath, 'rt')
            else:
                file_handle = open(filepath, 'rt')

            for line in file_handle:
                if line.startswith('#') or line.startswith('track'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < max(4, id_column + 1):
                    continue

                chr_name = fields[self.BED_INDEX['chr']]
                # BED format is 0-based, convert to 1-based
                start = int(fields[self.BED_INDEX['coord_start']]) + 1
                end = int(fields[self.BED_INDEX['coord_end']])

                # Generate feature ID based on feature type
                if feature_type == 'promoter':
                    # EPD promoters use coordinate-based IDs: SO:{chr_start_end_assembly}
                    region_id = build_regulatory_region_id(chr_name, start, end)
                    feature_id = f"SO:{region_id}"
                elif feature_type == 'non_coding_rna':
                    # RNAcentral format: URS000035F234_9606 -> RNACENTRAL:URS000035F234
                    if id_column < len(fields):
                        raw_id = fields[id_column].strip()
                        if raw_id and '_' in raw_id:
                            feature_id = f"RNACENTRAL:{raw_id.split('_')[0]}"
                        else:
                            feature_id = raw_id
                    else:
                        continue
                else:
                    # Generic BED feature
                    if id_column < len(fields):
                        feature_id = fields[id_column].strip()
                    else:
                        continue

                if feature_id:
                    yield (feature_id, chr_name, start, end)

            file_handle.close()

        except Exception as e:
            print(f"Error processing BED file for {feature_type}: {e}")

    def _process_promoter_features(self, filepath):
        """
        Process promoter features from ENCODE cCRE format.

        Yields:
            Tuples of (feature_id, chr, start, end)
        """
        try:
            if filepath.endswith('.gz'):
                file_handle = gzip.open(filepath, 'rt')
            else:
                file_handle = open(filepath, 'rt')

            for line in file_handle:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue

                element_type = fields[5]
                # Check if this is a promoter element
                is_promoter = element_type.startswith("PLS") or "pls" in element_type.lower()

                if not is_promoter:
                    continue

                chr_name = fields[0]
                # BED format is 0-based, convert to 1-based
                start = int(fields[1]) + 1
                end = int(fields[2]) + 1

                # Create promoter ID using SO:0000167 (promoter) prefix
                feature_id = f"SO:0000167:{chr_name}_{start}_{end}"

                yield (feature_id, chr_name, start, end)

            file_handle.close()

        except Exception as e:
            print(f"Error processing promoter file: {e}")

    def get_nodes(self):
        """This adapter only creates edges, not nodes."""
        return
        yield  # Make this a generator

    def get_edges(self):
        """
        Generate edges between genomic features and structural variants.

        Yields:
            Tuples of (source_id, target_id, label, properties)
        """
        edge_count = 0

        # Process each feature type
        for feature_type, source_config in self.feature_sources.items():
            if not source_config or 'filepath' not in source_config:
                continue

            filepath = source_config['filepath']
            file_format = source_config.get('format', 'gtf')

            print(f"Processing {feature_type} features from {filepath}")

            # Get feature iterator based on format
            if file_format == 'gtf':
                feature_iterator = self._process_gtf_features(filepath, feature_type)
            elif file_format == 'bed':
                id_column = source_config.get('id_column', 3)
                feature_iterator = self._process_bed_features(filepath, feature_type, id_column)
            elif file_format == 'promoter':
                feature_iterator = self._process_promoter_features(filepath)
            else:
                print(f"Unsupported format '{file_format}' for {feature_type}")
                continue

            # Process each feature
            for feature_id, chr_name, start, end in feature_iterator:
                # Get structural variants on the same chromosome
                sv_list = self.structural_variants.get(chr_name, [])

                # Compare with each structural variant
                for sv in sv_list:
                    sv_id = sv['id']
                    sv_start = sv['start']
                    sv_end = sv['end']

                    # Check relationship type based on label
                    should_create_edge = False

                    if self.label == 'located_in':
                        # Feature is fully contained within SV
                        should_create_edge = self._is_fully_contained(
                            chr_name, start, end,
                            chr_name, sv_start, sv_end
                        )

                    elif self.label == 'overlaps_with':
                        # Feature overlaps but is not fully contained
                        should_create_edge = self._has_overlap_but_not_contained(
                            chr_name, start, end,
                            chr_name, sv_start, sv_end
                        )

                    if should_create_edge:
                        props = {}

                        if self.write_properties:
                            props['chr'] = chr_name
                            props['feature_start'] = start
                            props['feature_end'] = end
                            props['sv_start'] = sv_start
                            props['sv_end'] = sv_end

                            if self.add_provenance:
                                props['source'] = self.source
                                props['sv_source'] = sv['source']

                        # Create edge from feature to structural variant
                        yield feature_id, sv_id, self.label, props
                        edge_count += 1

        print(f"Generated {edge_count} {self.label} edges")
