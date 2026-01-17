"""Smoke test for CORAL package.

This test verifies that the package can be imported and basic functionality works.
It does not download large datasets or run full pipelines.
"""

import os
import sys
import tempfile
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def test_imports():
    """Test that main modules can be imported."""
    from coral import MutationExtractionPipeline, MultiSpeciesMutationPipeline
    from coral import utils, cleanup_manager, genome_manager, alignment_manager
    from coral import pileup_manager, mutation_extractor_manager
    from coral import multiple_species_mutation_extractor_manager, plot_utils
    from coral import run_phylip, multiple_species_utils
    
    assert MutationExtractionPipeline is not None
    assert MultiSpeciesMutationPipeline is not None

def test_cli_import():
    """Test that CLI can be imported."""
    from coral.cli import main
    assert callable(main)

def test_utils():
    """Test utility functions."""
    from coral.utils import log, get_top_n_chromosomes
    
    # Test log function (should not crash)
    log("Test message", verbose=True)
    log("Test message", verbose=False)
    
    # Test get_top_n_chromosomes with a temporary fai file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fai', delete=False) as f:
        f.write("chr1\t1000\t0\t60\t61\n")
        f.write("chr2\t2000\t1000\t60\t61\n")
        f.write("chr3\t500\t3000\t60\t61\n")
        fai_path = f.name
    
    try:
        top_chroms = get_top_n_chromosomes(fai_path, n=2)
        assert len(top_chroms) == 2
        assert top_chroms[0] == "chr2"  # Longest
        assert top_chroms[1] == "chr1"
    finally:
        os.unlink(fai_path)

def test_pipeline_initialization():
    """Test that pipeline classes can be initialized."""
    from coral import MutationExtractionPipeline, MultiSpeciesMutationPipeline
    
    # Test single pipeline initialization
    pipeline = MutationExtractionPipeline(
        species_list=[("Species1", "ACC1"), ("Species2", "ACC2")],
        outgroup=("Outgroup", "OUT_ACC"),
        base_output_dir="./test_output"
    )
    assert pipeline.run_id is not None
    assert pipeline.output_dir is not None
    
    # Test multi-species pipeline initialization
    multi_pipeline = MultiSpeciesMutationPipeline(
        species_list=[("Species1", "ACC1"), ("Species2", "ACC2"), ("Outgroup", "OUT_ACC")],
        base_output_dir="./test_output",
        outgroup="Outgroup"
    )
    assert multi_pipeline.run_id is not None
    assert multi_pipeline.output_dir is not None

if __name__ == "__main__":
    test_imports()
    test_cli_import()
    test_utils()
    test_pipeline_initialization()
    print("All smoke tests passed!")

