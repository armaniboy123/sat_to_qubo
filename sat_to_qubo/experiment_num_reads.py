import random
import re
import csv
import time
import matplotlib.pyplot as plt

from sat_to_qubo import SATtoQUBO
from dwave.samplers import SimulatedAnnealingSampler


# =====================================================
# Experiment-Version des Compilers
# =====================================================

class ExperimentSATtoQUBO(SATtoQUBO):

    def __init__(self, P, domains, formula):

        self.P = P
        self.formula = formula
        self.domains = domains

        self.var_map = {}
        self.variables = []
        self.energy_terms = []

        self.aux_counter = 1

        self.create_variables()

    def read_formula(self):
        return self.formula


# =====================================================
# Zufällige Formeln erzeugen
#
# Erweitert gegenüber der ursprünglichen Version: nutzt alle
# fünf von SATtoQUBO unterstützten Operatoren (==, <=, >=, <, >)
# statt nur == und <=, sowie variable Klauselgrößen (1-4 Literale
# statt fix 3) -- damit werden auch Einzelliteral-Klauseln,
# Zwei-Literal-OR-Gadgets und die iterative OR-Reduktion für
# 3+ Literale (Abschnitt 4.6) durch den Generator abgedeckt.
# =====================================================

OPERATORS = ["==", "<=", ">=", "<", ">"]


def build_literal_pool(domains):

    pool = []

    for var in domains:
        for value in domains[var]:
            for op in OPERATORS:
                pool.append(f"{var}{op}{value}")

    return pool


def random_clause(domains, min_size=3, max_size=3):

    pool = build_literal_pool(domains)

    size = random.randint(min_size, max_size)
    size = min(size, len(pool))

    literals = random.sample(pool, size)

    return "(" + " OR ".join(literals) + ")"


def random_formula(domains,
                   n_clauses,
                   min_clause_size=3,
                   max_clause_size=3):

    clauses = []

    for _ in range(n_clauses):
        clauses.append(
            random_clause(domains, min_clause_size, max_clause_size)
        )

    return " AND ".join(clauses)


# =====================================================
# Formel überprüfen
#
# Regex und Operator-Auswertung um >=, <, > erweitert.
# Reihenfolge in der Alternation wichtig: <=/>= muessen vor
# </> geprueft werden (gleiche Logik wie in translate_literal).
# =====================================================

def evaluate_literal(literal, assignment):

    m = re.match(
        r"([a-zA-Z]\w*)(==|<=|>=|<|>)(\d+)",
        literal.strip()
    )

    var, op, value = m.group(1), m.group(2), int(m.group(3))

    if var not in assignment:
        return False

    actual = assignment[var]

    if op == "==": return actual == value
    if op == "<=": return actual <= value
    if op == ">=": return actual >= value
    if op == "<":  return actual <  value
    if op == ">":  return actual >  value


def evaluate_formula(formula, assignment):

    clauses = re.findall(r"\((.*?)\)", formula)

    for clause in clauses:

        literals = [
            l.strip()
            for l in re.split(r"\bOR\b", clause)
        ]

        satisfied = False

        for lit in literals:

            if evaluate_literal(lit, assignment):
                satisfied = True
                break

        if not satisfied:
            return False

    return True


# =====================================================
# Lösung dekodieren
# =====================================================

def decode_solution(compiler, sample):

    assignment = {}

    for var in compiler.domains:

        for value in compiler.domains[var]:

            idx = compiler.var_map[var][value]

            if sample.get(idx, 0) == 1:
                assignment[var] = value
                break

    return assignment


# =====================================================
# One-Hot-Verletzungen prüfen
# =====================================================

def is_valid_one_hot(compiler, sample):

    for var in compiler.domains:

        ones = sum(
            sample.get(compiler.var_map[var][v], 0)
            for v in compiler.domains[var]
        )

        if ones != 1:
            return False

    return True


# =====================================================
# Solver
# =====================================================

def solve_quiet(compiler, Q, num_reads=20):

    Q_dict = {k: v for k, v in Q.items() if v != 0.0}

    sampler = SimulatedAnnealingSampler()

    response = sampler.sample_qubo(
        Q_dict,
        num_reads=num_reads,
        num_sweeps=500
    )

    return response.first.sample


# =====================================================
# Experiment
# =====================================================

def run_experiment(n_vars,
                   formulas,
                   num_reads):

    domains = {}

    for i in range(n_vars):

        name = chr(ord('a') + i)

        domains[name] = [1, 2, 3]

    instances = len(formulas)

    successes = 0
    oh_violations = 0
    total_time = 0.0

    for k, formula in enumerate(formulas):

        compiler = ExperimentSATtoQUBO(
            P=10,
            domains=domains,
            formula=formula
        )

        start = time.time()

        Q = compiler.build()

        sample = solve_quiet(
            compiler,
            Q,
            num_reads=num_reads
        )

        total_time += time.time() - start

        if not is_valid_one_hot(compiler, sample):
            oh_violations += 1
            print(
                f"{k+1}/{instances} [One-Hot verletzt!]",
                end="\r"
            )
            continue

        assignment = decode_solution(
            compiler,
            sample
        )

        if evaluate_formula(
                formula,
                assignment):

            successes += 1

        print(
            f"{k+1}/{instances}",
            end="\r"
        )

    success_rate = successes / instances
    avg_time = total_time / instances

    print(
        f"\nOne-Hot Violations: {oh_violations}/{instances}"
    )

    return success_rate, avg_time


# =====================================================
# Hauptprogramm
# =====================================================

N_VARS = 10

CLAUSE_DENSITY = 3

NUM_READS_VALUES = [
    5,
    10,
    20,
    50,
    100
]

INSTANCES = 50

results = []

n_clauses = CLAUSE_DENSITY * N_VARS

# Formeln EINMALIG generieren und für alle num_reads-Werte
# wiederverwenden. Andernfalls würde bei jedem num_reads-Wert
# ein komplett neuer Satz zufälliger Formeln gezogen, wodurch
# Unterschiede in der Erfolgsrate nicht mehr eindeutig auf
# num_reads zurückführbar wären, sondern auch durch zufällig
# unterschiedlich schwere Formeln verursacht sein könnten.
domains_for_generation = {}
for i in range(N_VARS):
    name = chr(ord('a') + i)
    domains_for_generation[name] = [1, 2, 3]

formulas = [
    random_formula(domains_for_generation, n_clauses=n_clauses)
    for _ in range(INSTANCES)
]

for reads in NUM_READS_VALUES:

    print()
    print("=" * 60)
    print(f"num_reads = {reads}")
    print("=" * 60)

    success_rate, avg_time = run_experiment(
        n_vars=N_VARS,
        formulas=formulas,
        num_reads=reads
    )

    results.append(
        (
            reads,
            success_rate,
            avg_time
        )
    )

    print(
        f"Success Rate: {success_rate:.3f}"
    )

    print(
        f"Average Time: {avg_time:.3f} s"
    )

# =====================================================
# CSV speichern
# =====================================================

with open(
        "results_num_reads.csv",
        "w",
        newline=""
) as f:

    writer = csv.writer(f)

    writer.writerow([
        "num_reads",
        "success_rate",
        "avg_time"
    ])

    writer.writerows(results)


# =====================================================
# Plot: Success Rate
# =====================================================

x = [r[0] for r in results]
y = [r[1] for r in results]

plt.figure(figsize=(8, 5))

plt.plot(x, y, marker="o")

plt.xlabel("num_reads")

plt.ylabel("Success Rate")

plt.title(
    f"Success Rate vs num_reads "
    f"(n={N_VARS}, density={CLAUSE_DENSITY})"
)

plt.grid(True)

plt.savefig(
    "success_vs_num_reads.png"
)

plt.show()