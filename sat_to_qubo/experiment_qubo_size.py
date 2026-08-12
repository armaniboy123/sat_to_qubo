import random
import matplotlib.pyplot as plt
import numpy as np 

from sat_to_qubo import SATtoQUBO


# =====================================================
# Experiment-Version
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
# Zufällige Formel
#
# Erweitert gegenüber der ursprünglichen Version: nutzt alle
# fünf von SATtoQUBO unterstützten Operatoren (==, <=, >=, <, >)
# statt nur == und <=, konsistent mit den Generatoren der
# anderen Experimente in diesem Kapitel. Die Klauselgröße
# bleibt fix bei 3 Literalen.
# =====================================================

OPERATORS = ["==", "<=", ">=", "<", ">"]


def build_literal_pool(domains):

    pool = []

    for var in domains:
        for value in domains[var]:
            for op in OPERATORS:
                pool.append(f"{var}{op}{value}")

    return pool


def random_clause(domains, clause_size=3):

    pool = build_literal_pool(domains)

    clause = random.sample(pool, clause_size)

    return "(" + " OR ".join(clause) + ")"


def random_formula(domains,
                   n_clauses,
                   clause_size=3):

    clauses = []

    for _ in range(n_clauses):
        clauses.append(
            random_clause(domains, clause_size)
        )

    return " AND ".join(clauses)


# =====================================================
# Experiment
# =====================================================

def run_experiment(n_vars,
                   density=3,
                   instances=30):

    domains = {}

    for i in range(n_vars):
        domains[f"x{i}"] = [1, 2, 3]

    n_clauses = density * n_vars

    qubo_sizes = []
    helper_sizes = []

    for _ in range(instances):

        formula = random_formula(
            domains,
            n_clauses
        )

        compiler = ExperimentSATtoQUBO(
            P=10,
            domains=domains,
            formula=formula
        )

        compiler.build()

        one_hot = sum(
            len(v)
            for v in compiler.domains.values()
        )

        helpers = compiler.n - one_hot

        helper_sizes.append(helpers)
        qubo_sizes.append(compiler.n)

    avg_qubo = sum(qubo_sizes) / instances
    avg_helper = sum(helper_sizes) / instances

    return avg_qubo, avg_helper


# =====================================================
# Hauptprogramm
# =====================================================

VARIABLES = [
    5,
    10,
    15,
    20,
    25,
    30
]

CLAUSE_DENSITY = 3

qubo = []
helpers = []

for n in VARIABLES:

    print(f"{n} Variablen...")

    q, h = run_experiment(
        n,
        density=CLAUSE_DENSITY,
        instances=30
    )

    qubo.append(q)
    helpers.append(h)


# =====================================================
# Durchschnittswerte ausgeben
# =====================================================

print()
print("=" * 70)
print("Durchschnittswerte")
print("=" * 70)

print(f"{'SAT Vars':>10} {'QUBO Vars':>12} {'Helper Vars':>14}")

for n, q, h in zip(VARIABLES, qubo, helpers):
    print(f"{n:>10} {q:>12.1f} {h:>14.1f}")


# =====================================================
# Lineare Näherung berechnen
# =====================================================

qubo_m, qubo_b = np.polyfit(VARIABLES, qubo, 1)
helper_m, helper_b = np.polyfit(VARIABLES, helpers, 1)

print()
print("=" * 70)
print("Geschätztes Verhältnis")
print("=" * 70)

print(f"QUBO(n)   ≈ {qubo_m:.2f} · n + {qubo_b:.2f}")
print(f"Helper(n) ≈ {helper_m:.2f} · n + {helper_b:.2f}")

print()
print("mit")
print("n          = Anzahl der SAT-Variablen")
print("QUBO(n)    = Anzahl der QUBO-Variablen")
print("Helper(n)  = Anzahl der Hilfsvariablen")

print()
print(f"⇒ Pro zusätzlicher SAT-Variable entstehen")
print(f"  durchschnittlich etwa {qubo_m:.2f} QUBO-Variablen.")


# =====================================================
# Plot Gesamtgröße
# =====================================================

plt.figure(figsize=(8, 5))

plt.plot(
    VARIABLES,
    qubo,
    marker="o",
    label="QUBO variables"
)

plt.xlabel("SAT variables")
plt.ylabel("Number of QUBO variables")
plt.title("Growth of QUBO size")

plt.grid(True)
plt.legend()

plt.savefig("qubo_growth.png")

plt.show()


# =====================================================
# Plot Hilfsvariablen
# =====================================================

plt.figure(figsize=(8, 5))

plt.plot(
    VARIABLES,
    helpers,
    marker="o",
    color="red",
    label="Helper variables"
)

plt.xlabel("SAT variables")
plt.ylabel("Helper variables")
plt.title("Growth of helper variables")

plt.grid(True)
plt.legend()

plt.savefig("helper_growth.png")

plt.show()