"""Tests for model.py. Each test checks one claim made in the README.
Run with:  pytest -v
"""
import pytest
import model as m

# Per-attractor basin counts from running the authors' rules file through BoolNet
BOOLNET_COUNTS = [551768, 646816, 54336, 72448, 200808, 91136,
                  107096, 107104, 36864, 36864, 4096, 4096,
                  3496, 24576, 24576, 53760, 65536, 11776]

MUTATIONS = {"foxo3 LOF": {"foxo3": 0}, "antiox GOF": {"antiox": 1},
             "mef2c LOF": {"mef2c": 0}, "akt GOF": {"akt": 1},
             "meis1 GOF": {"meis1": 1}, "runx1 GOF": {"runx1": 1}}


# ---------- Shared results: each condition is simulated once for all tests ----------
@pytest.fixture(scope="module")
def wild_type():
    return m.attractor_data({})

@pytest.fixture(scope="module")
def mutants():
    return {label: m.attractor_data(fixed) for label, fixed in MUTATIONS.items()}

def pct_included(counts, cell):
    """Percentage out of all states, cycles included."""
    return 100 * counts.get(cell, 0) / sum(counts.values())

def pct_excluded(counts, cell):
    """Percentage out of states that reach a cell type, cycles excluded."""
    return 100 * counts.get(cell, 0) / (sum(counts.values()) - counts["NS"])


# ---------- Validation 1: the published cell states ----------
def test_eleven_figure_states_are_fixed_points():
    unchanged = [name for name, s in m.figure_states.items() if m.update_state(s) == s]
    assert len(unchanged) == 11

def test_gmp1_and_gmp3_differ_only_at_runx1():
    for name in ["gmp1", "gmp3"]:
        state = m.figure_states[name]
        new_state = m.update_state(state)
        changed = [n for n in m.NODES if new_state[n] != state[n]]
        assert changed == ["runx1"], name


# ---------- Validation 2: fixed points from scratch ----------
def test_exactly_13_fixed_points():
    assert len(m.find_fixed_points()) == 13


# ---------- Validation 3: basin sizes ----------
def test_wild_type_has_18_attractors(wild_type):
    assert len(wild_type) == 18

def test_wild_type_counts_every_state(wild_type):
    assert sum(wild_type.values()) == 2 ** 21

def test_wild_type_matches_boolnet(wild_type):
    assert sorted(wild_type.values()) == sorted(BOOLNET_COUNTS)

def test_wild_type_percentages_match_paper(wild_type):
    counts = m.basin_count(wild_type)
    expected = {"HSC": 0.18, "MEP": 84.37, "GMP": 4.27, "CLP": 11.17}
    for cell, pct in expected.items():
        assert round(pct_excluded(counts, cell), 2) == pct, cell


# ---------- Mutants: results reported in the paper ----------
def test_mutant_attractor_counts(mutants):
    expected = {"foxo3 LOF": 18, "antiox GOF": 18, "mef2c LOF": 14,
                "akt GOF": 16, "meis1 GOF": 18, "runx1 GOF": 16}
    for label, n in expected.items():
        assert len(mutants[label]) == n, label

def test_mutant_percentages_match_paper(mutants):
    checks = [("foxo3 LOF", "HSC", 0.067),
              ("antiox GOF", "HSC", 0.331),
              ("mef2c LOF", "CLP", 6.134)]
    for label, cell, pct in checks:
        counts = m.basin_count(mutants[label])
        assert round(pct_included(counts, cell), 3) == pct, label

def test_akt_gof_has_no_hsc(mutants):
    assert m.basin_count(mutants["akt GOF"]).get("HSC", 0) == 0

def test_table2_reductions_under_papers_method(wild_type, mutants):
    wt = m.basin_count(wild_type)
    for label, cell in [("meis1 GOF", "MEP"), ("runx1 GOF", "MEP"), ("runx1 GOF", "HSC")]:
        counts = m.basin_count(mutants[label])
        assert pct_included(counts, cell) < pct_excluded(wt, cell), (label, cell)


# ---------- The comparison findings ----------
def test_paper_fold_change_is_consistent_times_wt_fate_fraction(wild_type, mutants):
    wt = m.basin_count(wild_type)
    wt_fate_fraction = (sum(wt.values()) - wt["NS"]) / sum(wt.values())
    for label in mutants:
        counts = m.basin_count(mutants[label])
        for cell in ["HSC", "MEP", "GMP", "CLP"]:
            mut = pct_included(counts, cell)
            if mut == 0:
                continue
            paper_fc = mut / pct_excluded(wt, cell)
            consistent_fc = mut / pct_included(wt, cell)
            assert paper_fc / consistent_fc == pytest.approx(wt_fate_fraction), (label, cell)

def test_gmp_unchanged_under_consistent_method(wild_type, mutants):
    wt = m.basin_count(wild_type)
    for label in ["foxo3 LOF", "antiox GOF", "akt GOF"]:
        counts = m.basin_count(mutants[label])
        assert pct_included(counts, "GMP") == pytest.approx(pct_included(wt, "GMP")), label

def test_cycles_excluded_reverses_two_claims(wild_type, mutants):
    wt = m.basin_count(wild_type)
    meis1 = m.basin_count(mutants["meis1 GOF"])
    runx1 = m.basin_count(mutants["runx1 GOF"])
    assert pct_excluded(meis1, "MEP") / pct_excluded(wt, "MEP") == pytest.approx(1.005, abs=1e-3)
    assert pct_excluded(runx1, "HSC") / pct_excluded(wt, "HSC") == pytest.approx(1.873, abs=1e-3)

def test_akt_gof_total_states(mutants):
    assert sum(mutants["akt GOF"].values()) == 1048576