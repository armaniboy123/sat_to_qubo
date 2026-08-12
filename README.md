# SAT-to-QUBO

A Python framework for translating Boolean satisfiability problems over **categorical, finite-domain variables** into **Quadratic Unconstrained Binary Optimization (QUBO)** formulations, suitable for quantum annealing hardware or classical simulated annealing.

## Repository Structure

```
.
├── README.md
└── sat_to_qubo/
    ├── sat_to_qubo.py               # core translation framework
    ├── experiment_clause_density.py # Influence of Clause Density
    ├── experiment_variables.py      # Influence of the Number of SAT Variables
    ├── experiment_num_reads.py      # Influence of num_reads
    └── experiment_qubo_size.py      # Growth of the QUBO Representation
```

All commands below assume you are running them from inside the `sat_to_qubo/` directory:

```bash
cd sat_to_qubo
```

## Requirements

- Python 3.9 or later
- [`dwave-samplers`](https://pypi.org/project/dwave-samplers/) (provides `SimulatedAnnealingSampler`)
- `numpy`
- `matplotlib` (only required for the experiment scripts)

### Installation

```bash
pip install dwave-samplers numpy matplotlib
```

---

## 1. Running the Core Framework (`sat_to_qubo.py`)

The core framework can be run directly and interactively from the command line:

```bash
python3 sat_to_qubo.py
```

You will be prompted for three inputs, in order:

### Step 1 — Penalty value

```
Penalty:
> 10
```

Enter a positive integer. `P = 10` is a reasonable default for small examples.

### Step 2 — Variable domains

```
Wertebereiche:
Format: a=1-5 b=1-5 c=1-5
> a=1-4 b=1-4 c=1-4
```

Declare every categorical variable you intend to use in the formula, along with its finite domain, as `name=low-high`. All variables used in the formula **must** be declared here first, or the parser will raise an error.

### Step 3 — The SAT formula

```
SAT Formel:
Format: (a==3 OR b<=2) AND (c==1 OR d==2 OR a==4)
Operatoren: ==, <=, >=, <, >
> (a==2 OR b<=1) AND (c>=3)
```

The formula must be given in conjunctive normal form: a conjunction (`AND`) of clauses, where each clause is a disjunction (`OR`) of literals. Each literal has the form `variable<operator>value`, using one of the five supported operators: `==`, `<=`, `>=`, `<`, `>`. Clauses may contain any number of literals (including just one).

### Output

The program prints, in order:

1. **Energy function** — a human-readable listing of every penalty term added to the objective.
2. **QUBO matrix** — the full matrix, labeled with variable names (One-Hot variables and auxiliary variables `h`/`z`/`u`). In addition to printing to the console, the labeled matrix is written unabridged (without terminal line wrapping) to `output.txt` in the current directory — useful for instances with many variables, where the console output wraps and becomes hard to read.
3. **Variable list** — mapping from index to variable name.
4. **Solver output** — the best solution found by simulated annealing, including:
   - the raw binary assignment,
   - the decoded categorical assignment (with a warning if any One-Hot constraint was violated),
   - the original formula with the found assignment substituted in (e.g. `(TRUE OR FALSE) AND (TRUE)`), together with the overall result (`TRUE`/`FALSE`).

### Using `SATtoQUBO` programmatically

The class can also be used without the interactive prompts, e.g. from a script or notebook:

```python
from sat_to_qubo import SATtoQUBO

compiler = SATtoQUBO.__new__(SATtoQUBO)
compiler.P = 10
compiler.domains = {"a": [1, 2, 3], "b": [1, 2, 3]}
compiler.var_map = {}
compiler.variables = []
compiler.energy_terms = []
compiler.aux_counter = 1
compiler.create_variables()

formula = "(a==2 OR b<=1)"
compiler.formula = formula
compiler.add_one_hot()
compiler.parse_formula(formula)

response = compiler.solve(compiler.Q, num_reads=100)
```

---

## 2. Running the Experiments

Each experiment script is self-contained and can be run directly:

```bash
python3 experiment_clause_density.py
python3 experiment_variables.py
python3 experiment_num_reads.py
python3 experiment_qubo_size.py
```

All four scripts import the core framework from `sat_to_qubo.py`, so this file must be located in the same directory.

Each script prints progress to the console, writes its results to a `.csv` file, and (except `experiment_qubo_size.py`, which prints a summary table instead) saves one or more `.png` plots — all in the current working directory.

### 2.1 `experiment_clause_density.py`

Measures success rate and average runtime as the clause density (`clauses / variable`) is varied, at a fixed number of variables.

| Parameter | Default | Description |
|---|---|---|
| `N_VARS` | `10` | Number of SAT variables (fixed) |
| `CLAUSE_DENSITIES` | `[1, 2, 3, 4, 5, 6]` | Clause-to-variable ratios tested |
| `INSTANCES` | `50` | Random formulas generated per configuration |

**Output:** `results_clause_density.csv`, `success_vs_density.png`, `runtime_vs_density.png`

### 2.2 `experiment_variables.py`

Measures success rate and average runtime as the number of SAT variables is varied, at a fixed clause density.

| Parameter | Default | Description |
|---|---|---|
| `VARIABLES` | `[5, 10, 15, 20, 25]` | Numbers of variables tested |
| `CLAUSE_DENSITY` | `2` | Clauses per variable (fixed) |
| `INSTANCES` | `50` | Random formulas generated per configuration |

**Output:** `results_variables.csv`, `success_vs_variables.png`, `runtime_vs_variables.png`

### 2.3 `experiment_num_reads.py`

Measures how the solver's success rate depends on the number of independent reads (`num_reads`), at fixed variable count and clause density. To isolate the effect of `num_reads`, the **same** set of random formulas is reused across all tested values.

| Parameter | Default | Description |
|---|---|---|
| `N_VARS` | `10` | Number of SAT variables (fixed) |
| `CLAUSE_DENSITY` | `3` | Clauses per variable (fixed) |
| `NUM_READS_VALUES` | `[5, 10, 20, 50, 100]` | Values of `num_reads` tested |
| `INSTANCES` | `50` | Random formulas generated once and reused for every `num_reads` value |

**Output:** `results_num_reads.csv`, `success_vs_num_reads.png`

### 2.4 `experiment_qubo_size.py`

Measures the number of QUBO variables (total and auxiliary) as a function of the number of SAT variables, without solving the generated instances. Also fits and prints a linear regression.

| Parameter | Default | Description |
|---|---|---|
| `VARIABLES` | `[5, 10, 15, 20, 25, 30]` | Numbers of variables tested |
| `CLAUSE_DENSITY` | `3` | Clauses per variable (fixed) |
| `instances` (in `run_experiment`) | `30` | Random formulas averaged per configuration |

**Output:** printed summary table and regression coefficients (`qubo_growth.png`, `helper_growth.png`)

---

## 3. Adjusting Experiment Parameters

All experiment parameters (variable ranges, clause density, penalty value `P`, `num_sweeps`, `num_reads`, number of instances) are defined as constants near the top of the "Hauptprogramm" (main program) section of each script, and can be edited directly to reproduce results under different configurations.

---

## 4. Notes

- All experiments use `SimulatedAnnealingSampler` from `dwave-samplers` as a classical proxy for quantum annealing hardware; no D-Wave account or hardware access is required.
- The random formula generator used in the experiment scripts draws literals using all five supported operators (`==`, `<=`, `>=`, `<`, `>`) and, unless stated otherwise above, uses a fixed clause size of three literals.
- Running an experiment script will overwrite any previously generated `.csv` and `.png` files of the same name in the working directory.
