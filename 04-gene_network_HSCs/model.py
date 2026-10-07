import itertools, numpy as np, hashlib, inspect, os, json, pandas as pd

#Genere gulatory  network with asynchronous updated produced patterns (Figure 3: Herrera et al)
hsc = {"oxygen":0, "runx1":0, "meis1":1, "hif1":1, "foxo3":1, "p53":1, "gata2":1, "gata1":0, "pu1":0, "cebpa":0, "ikzf1":1, "gfi1":0, "mef2c":0, "mtor":0, "ampk":1, "akt":0, "h2o2":0, "o2":0, "sod":1, "antiox":1, "oxphos":0}
mep1 = {"oxygen":0, "runx1":0, "meis1":0, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
mep2 = {"oxygen":0, "runx1":0, "meis1":1, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
mep3 = {"oxygen":0, "runx1":1, "meis1":1, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
mep4 = {"oxygen":1, "runx1":1, "meis1":1, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}
mep5 = {"oxygen":1, "runx1":0, "meis1":1, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}
mep6 = {"oxygen":1, "runx1":0, "meis1":0, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":1, "pu1":0, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}
gmp1 = {"oxygen":0, "runx1":0, "meis1":1, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":1, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":1, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
gmp2 = {"oxygen":0, "runx1":1, "meis1":1, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":1, "cebpa":1, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
gmp3 = {"oxygen":1, "runx1":0, "meis1":1, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":1, "cebpa":0, "ikzf1":0, "gfi1":0, "mef2c":1, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}
gmp4 = {"oxygen":1, "runx1":1, "meis1":1, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":1, "cebpa":1, "ikzf1":0, "gfi1":0, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}
clp1 = {"oxygen":0, "runx1":0, "meis1":0, "hif1":1, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":0, "cebpa":0, "ikzf1":1, "gfi1":1, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":0}
clp2 = {"oxygen":1, "runx1":0, "meis1":0, "hif1":0, "foxo3":0, "p53":0, "gata2":0, "gata1":0, "pu1":0, "cebpa":0, "ikzf1":1, "gfi1":1, "mef2c":0, "mtor":1, "ampk":0, "akt":1, "h2o2":1, "o2":1, "sod":0, "antiox":0, "oxphos":1}

figure_states = {"hsc": hsc, "mep1": mep1, "mep2": mep2, "mep3": mep3,
                 "mep4": mep4, "mep5": mep5, "mep6": mep6,
                 "gmp1": gmp1, "gmp2": gmp2, "gmp3": gmp3, "gmp4": gmp4,
                 "clp1": clp1, "clp2": clp2}

NODES = ["oxygen", "runx1", "meis1", "hif1", "foxo3", "p53", "gata2",
         "gata1", "pu1", "cebpa", "ikzf1", "gfi1", "mef2c", "mtor",
         "ampk", "akt", "h2o2", "o2", "sod", "antiox", "oxphos"]

fixed = [{"foxo3": 0}, {"antiox": 1}, {"mef2c": 0}, {"akt": 1},
         {"meis1": 1}, {"runx1": 1}, {}]

CELL_TYPES = ["HSC", "MEP", "GMP", "CLP", "NS"]


#Function calculating every node's next value from the current state using the rules declared in the paper
def update_state(cs, fixed=None):
    ns = {}
    ns["oxygen"] = cs["oxygen"]
    ns["runx1"] = int((cs["pu1"] or cs["gata2"] or cs["runx1"]) and not cs["ikzf1"])
    ns["meis1"] = int((cs["meis1"] or cs["runx1"] or cs["pu1"]) and not cs["gfi1"])
    ns["hif1"] = int((not cs["oxygen"] and cs["o2"] and not (cs["p53"] or cs["foxo3"] or cs["ampk"])) or (not cs["oxygen"] and cs["meis1"]))
    ns["foxo3"] = int((cs["ampk"] and (cs["hif1"] or cs["p53"])) and not (cs["akt"] or cs["cebpa"]))
    ns["p53"] = int(cs["hif1"] and not (cs["akt"] or cs["gata1"] or cs["pu1"]))
    ns["gata2"] = int((cs["gata2"] or cs["p53"]) and not (cs["gata1"] or cs["pu1"] or cs["gfi1"]))
    ns["gata1"] = int((cs["gata1"] or cs["gata2"] or cs["akt"] or cs["runx1"]) and not (cs["pu1"] or cs["p53"] or cs["ikzf1"]))
    ns["pu1"] = int((cs["pu1"] or cs["runx1"] or (cs["cebpa"] and cs["ikzf1"])) and not (cs["gata1"] or cs["gata2"] or cs["gfi1"]))
    ns["cebpa"] = int((cs["pu1"] and cs["runx1"]) and not cs["mef2c"])
    ns["ikzf1"] = int((cs["mef2c"] or cs["ikzf1"] or cs["runx1"]) and not (cs["cebpa"] or cs["pu1"] or cs["gata1"]))
    ns["gfi1"] = int((cs["ikzf1"] or cs["cebpa"]) and not (cs["pu1"] or cs["p53"]))
    ns["mef2c"] = int(cs["pu1"] and not cs["cebpa"])
    ns["mtor"] = int(cs["akt"] and not cs["ampk"] and not cs["p53"])
    ns["ampk"] = int((not cs["akt"] or not cs["oxphos"] or not cs["hif1"]) and cs["p53"])
    ns["akt"] = int(cs["h2o2"] or not (cs["p53"] or cs["foxo3"]))
    ns["h2o2"] = int((cs["gata1"] or cs["hif1"] or cs["pu1"] or cs["oxphos"] or (cs["sod"] and cs["o2"])) and not cs["antiox"])
    ns["o2"] = int(not cs["sod"] or cs["oxphos"])
    ns["sod"] = int(cs["p53"] or cs["foxo3"])
    ns["antiox"] = int(cs["p53"] or cs["foxo3"])
    ns["oxphos"] = int((cs["mtor"] or cs["akt"]) and not (cs["foxo3"] or cs["hif1"]))

    if fixed:
        for gene, value in fixed.items():
            ns[gene] = value
    return ns

def check_figure_states():
    for name, state in figure_states.items():
        new_state = update_state(state)
        if new_state == state:
            print(f"{name}: unchanged (fixed point)")
        else:
            changed = [n for n in NODES if new_state[n] != state[n]]
            print(f"{name}: changed at {changed}")

def find_fixed_points():
    fixed_points = []
    for i, values in enumerate(itertools.product([0, 1], repeat=len(NODES))):
        state = dict(zip(NODES, values))

        
        if update_state(state) == state:
            fixed_points.append(state)

        
        if i % 500000 == 0:
            print(f"Checked {i:,} states, fixed points so far: {len(fixed_points)}")

    print(f"\nTotal fixed points found: {len(fixed_points)}")
    return fixed_points

def cell_state_check(fixed_points):
    for fp in fixed_points:
        exact = [name for name, s in figure_states.items() if s == fp]
        if exact:
            print(f"Exact match: {exact[0]}")
        else:
            # Find the Figure 3 column that differs in the fewest nodes
            name, closest = min(figure_states.items(),
                                key=lambda kv: sum(kv[1][n] != fp[n] for n in NODES))
            differences = [n for n in NODES if closest[n] != fp[n]]
            print(f"No exact match. Closest: {name}, differs at: {differences}")


def to_key(state):
    """Turn a state dictionary into a 0/1 string in NODES order (same format as BoolNet)."""
    return "".join(str(state[n]) for n in NODES)

def attractor_data(fixed):
    """For every starting state, follow the updates until a state repeats,
    then count which attractor it ended in."""
    attractors = {}
    for values in itertools.product([0, 1], repeat=len(NODES)):
        state = dict(zip(NODES, values))   # starting state

        # A mutated gene can't start at the wrong value, so skip those states
        if any(state[g] != v for g, v in fixed.items()):
            continue

        temp_states = []
        while True:
            temp_states.append(state)
            state = update_state(state, fixed)

            if state in temp_states:                            # a state repeated: we've reached an attractor
                first = temp_states.index(state)                # where the repeat first appeared
                attractor_states = temp_states[first:]          # the attractor: repeat point to the end
                name = "|".join(sorted(to_key(s) for s in attractor_states))
                attractors[name] = attractors.get(name, 0) + 1  # count this starting state
                break
    return attractors

def cell_type(name):
    """Label an attractor by its marker gene. Cycles (names containing '|') get 'NS'."""
    if "|" in name:
        return "NS"

    # Turn the 0/1 string back into a dictionary so we can check genes by name
    state = dict(zip(NODES, (int(c) for c in name)))

    if state["gata2"]:
        return "HSC"
    if state["gata1"]:
        return "MEP"
    if state["pu1"]:
        return "GMP"
    if state["gfi1"]:
        return "CLP"
    return "unknown"

def basin_count(attractors):
    type_counts = {}
    for name, count in attractors.items():
        label = cell_type(name)
        type_counts[label] = type_counts.get(label, 0) + count
    return type_counts

RULES_HASH = hashlib.md5(inspect.getsource(update_state).encode()).hexdigest()[:8]

def condition_label(fixed):
    if not fixed:
        return "wild_type"
    gene, value = next(iter(fixed.items()))
    return f"{gene}={value}"

def attractor_data_cached(fixed, folder="basin_cache"):
    os.makedirs(folder, exist_ok=True)
    fname = f"{condition_label(fixed).replace('=', '_')}_{RULES_HASH}.json"
    path = os.path.join(folder, fname)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    result = attractor_data(fixed)
    with open(path, "w") as f:
        json.dump(result, f)
    return result



def create_table(mutations):
    results = {}
    for fixed in mutations:
        label = condition_label(fixed)
        counts = basin_count(attractor_data_cached(fixed))
        assert "unknown" not in counts, f"{label}: attractor with no marker gene"
        counts = {c: counts.get(c, 0) for c in CELL_TYPES}
        expected = 2 ** (len(NODES) - len(fixed))
        assert sum(counts.values()) == expected, f"{label}: counted {sum(counts.values())}, expected {expected}"
        results[label] = counts
    results = {"wild_type": results.pop("wild_type"), **results}   # wild type first

    wt = results["wild_type"]
    wt_total = sum(wt.values())
    wt_fate = wt_total - wt["NS"]
    wt_ci = {c: 100 * wt[c] / wt_total for c in CELL_TYPES}                   # cycles included
    wt_ce = {c: 100 * wt[c] / wt_fate for c in CELL_TYPES if c != "NS"}       # cycles excluded

    rows = []
    for label, counts in results.items():
        total = sum(counts.values())
        fate = total - counts["NS"]
        for c in CELL_TYPES:
            n = counts[c]
            pct_ci = 100 * n / total
            pct_ce = np.nan if c == "NS" else 100 * n / fate
            row = {"Condition": label, "Cell type": c, "Basin size": n,
                   "Percent (CI)": pct_ci, "Percent (CE)": pct_ce,
                   "FC paper (mutant CI / WT CE)": np.nan,
                   "FC CI/CI": np.nan, "FC CE/CE": np.nan}
            if label != "wild_type":
                row["FC CI/CI"] = pct_ci / wt_ci[c]            # NS row shows how much cycles grew
                if c != "NS":
                    row["FC paper (mutant CI / WT CE)"] = pct_ci / wt_ce[c]
                    row["FC CE/CE"] = pct_ce / wt_ce[c]
            rows.append(row)                                    # zero basins give 0.0, not NaN

    df = pd.DataFrame(rows)
    df[df.columns[3:]] = df[df.columns[3:]].round(3)
    return df

def check_half_cycles(df):
    """If a mutant has exactly half the wild-type cycle states, CI/CI must equal CE/CE."""
    ns = df[df["Cell type"] == "NS"].set_index("Condition")["Basin size"]
    for cond in ns.index.drop("wild_type"):
        if 2 * ns[cond] == ns["wild_type"]:
            sub = df[(df["Condition"] == cond) & (df["Cell type"] != "NS")]
            assert np.allclose(sub["FC CI/CI"], sub["FC CE/CE"], atol=2e-3), cond
            print(f"{cond}: CI/CI == CE/CE  OK")