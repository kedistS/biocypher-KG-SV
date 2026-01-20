#!/usr/bin/env python3
"""
Standalone test script for coordinate comparison logic.
Tests the core logic without requiring full adapter imports.
"""


def _is_fully_contained(feat_chr, feat_start, feat_end, sv_chr, sv_start, sv_end):
    """Check if a feature is fully contained within a structural variant."""
    if feat_chr != sv_chr:
        return False
    return feat_start >= sv_start and feat_end <= sv_end


def _has_overlap_but_not_contained(feat_chr, feat_start, feat_end, sv_chr, sv_start, sv_end):
    """Check if a feature overlaps with a structural variant but is not fully contained."""
    if feat_chr != sv_chr:
        return False
    has_overlap = not (feat_end < sv_start or feat_start > sv_end)
    if not has_overlap:
        return False
    is_contained = feat_start >= sv_start and feat_end <= sv_end
    return not is_contained


def test_coordinate_comparison():
    """Test the coordinate comparison logic."""
    print("Testing coordinate comparison logic...")

    # Test 1: Feature fully contained within SV
    assert _is_fully_contained(
        'chr1', 1000, 2000,  # Feature
        'chr1', 500, 3000    # SV
    ) == True, "Test 1 failed: Feature should be fully contained"

    # Test 2: Feature not contained (different chromosome)
    assert _is_fully_contained(
        'chr1', 1000, 2000,  # Feature
        'chr2', 500, 3000    # SV
    ) == False, "Test 2 failed: Different chromosomes should not match"

    # Test 3: Feature partially overlapping (extends beyond SV start)
    assert _is_fully_contained(
        'chr1', 400, 1500,   # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 3 failed: Partially overlapping feature should not be contained"

    # Test 4: Feature partially overlapping (extends beyond SV end)
    assert _is_fully_contained(
        'chr1', 2500, 3500,  # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 4 failed: Partially overlapping feature should not be contained"

    # Test 5: Exact match (feature = SV)
    assert _is_fully_contained(
        'chr1', 1000, 2000,  # Feature
        'chr1', 1000, 2000   # SV
    ) == True, "Test 5 failed: Exact match should be contained"

    print("✓ All _is_fully_contained tests passed")

    # Test overlap detection (but not contained)

    # Test 6: Feature overlaps but extends beyond SV start
    assert _has_overlap_but_not_contained(
        'chr1', 400, 1500,   # Feature
        'chr1', 500, 3000    # SV
    ) == True, "Test 6 failed: Should overlap but not be contained"

    # Test 7: Feature overlaps but extends beyond SV end
    assert _has_overlap_but_not_contained(
        'chr1', 2500, 3500,  # Feature
        'chr1', 500, 3000    # SV
    ) == True, "Test 7 failed: Should overlap but not be contained"

    # Test 8: SV fully contained within feature (reverse containment)
    assert _has_overlap_but_not_contained(
        'chr1', 100, 5000,   # Feature
        'chr1', 500, 3000    # SV
    ) == True, "Test 8 failed: SV contained in feature should count as overlap"

    # Test 9: Feature fully contained (should NOT trigger overlap_but_not_contained)
    assert _has_overlap_but_not_contained(
        'chr1', 1000, 2000,  # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 9 failed: Fully contained should not trigger overlap"

    # Test 10: No overlap (feature before SV)
    assert _has_overlap_but_not_contained(
        'chr1', 100, 400,    # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 10 failed: No overlap should return False"

    # Test 11: No overlap (feature after SV)
    assert _has_overlap_but_not_contained(
        'chr1', 3500, 4000,  # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 11 failed: No overlap should return False"

    # Test 12: Adjacent regions (no overlap)
    assert _has_overlap_but_not_contained(
        'chr1', 3001, 4000,  # Feature
        'chr1', 500, 3000    # SV
    ) == False, "Test 12 failed: Adjacent regions should not overlap"

    # Test 13: Different chromosomes (no overlap)
    assert _has_overlap_but_not_contained(
        'chr1', 1000, 2000,  # Feature
        'chr2', 500, 3000    # SV
    ) == False, "Test 13 failed: Different chromosomes should not overlap"

    print("✓ All _has_overlap_but_not_contained tests passed")


def test_edge_cases():
    """Test edge cases and boundary conditions."""
    print("\nTesting edge cases...")

    # Test single base pair feature
    assert _is_fully_contained(
        'chr1', 1500, 1500,  # Single base
        'chr1', 1000, 2000   # SV
    ) == True, "Single base contained test failed"

    # Test single base pair SV
    assert _is_fully_contained(
        'chr1', 1500, 1600,  # Feature
        'chr1', 1500, 1500   # Single base SV
    ) == False, "Single base SV test failed"

    # Test feature touching SV start
    assert _has_overlap_but_not_contained(
        'chr1', 500, 1500,   # Feature starts at SV start
        'chr1', 500, 3000    # SV
    ) == False, "Feature starting at SV boundary should be contained"

    # Test feature touching SV end
    assert _has_overlap_but_not_contained(
        'chr1', 2000, 3000,  # Feature ends at SV end
        'chr1', 500, 3000    # SV
    ) == False, "Feature ending at SV boundary should be contained"

    print("✓ All edge case tests passed")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Running Coordinate Comparison Logic Tests")
    print("=" * 60)

    tests_passed = 0
    tests_failed = 0

    try:
        test_coordinate_comparison()
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Coordinate comparison test failed: {e}")
        tests_failed += 1
    except Exception as e:
        print(f"✗ Unexpected error in coordinate comparison test: {e}")
        tests_failed += 1

    try:
        test_edge_cases()
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Edge case test failed: {e}")
        tests_failed += 1
    except Exception as e:
        print(f"✗ Unexpected error in edge case test: {e}")
        tests_failed += 1

    print("\n" + "=" * 60)
    print(f"Test Results: {tests_passed} passed, {tests_failed} failed")
    print("=" * 60)

    return tests_failed == 0


if __name__ == '__main__':
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)
