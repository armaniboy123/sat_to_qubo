import numpy as np
import re
from itertools import combinations
from dwave.samplers import SimulatedAnnealingSampler


class SATtoQUBO:

    def __init__(self):

        self.P = self.read_penalty()

        self.domains = {}
        self.var_map = {}
        self.variables = []

        self.energy_terms = []

        self.aux_counter = 1

    # ---------------------------------

    def read_penalty(self):

        print("\nPenalty:")

        return int(input("> "))

    # ---------------------------------

    def read_domains(self):

        print("\nWertebereiche:")
        print("Format: a=1-5 b=1-5 c=1-5")

        line = input("> ")

        for p in line.split():

            var, interval = p.split("=")

            low, high = map(
                int,
                interval.split("-")
            )

            self.domains[var] = list(
                range(low, high + 1)
            )

        self.create_variables()

    # ---------------------------------

    def read_formula(self):

        print("\nSAT Formel:")
        print("Format: (a==3 OR b<=2) AND (c==1 OR d==2 OR a==4)")
        print("Operatoren: ==, <=, >=, <, >")

        return input("> ")

    # ---------------------------------

    def create_variables(self):

        idx = 0

        for var in self.domains:

            self.var_map[var] = {}

            for val in self.domains[var]:

                self.var_map[var][val] = idx

                self.variables.append(
                    f"{var}_{val}"
                )

                idx += 1

        self.n = idx

        # Dynamisch wachsendes Dictionary statt fixer Matrixgröße.
        # Keys sind (i, j) mit i <= j; Werte sind die Q-Koeffizienten.
        self.Q = {}

    # ---------------------------------

    def create_aux(self, prefix="z"):

        name = f"{prefix}{self.aux_counter}"
        idx  = self.n

        self.variables.append(name)

        self.n           += 1
        self.aux_counter += 1

        return idx, name

    # ---------------------------------

    def add_linear(self, i, val):
        key = (i, i)
        self.Q[key] = self.Q.get(key, 0.0) + val

    def add_quad(self, i, j, val):
        if i == j:
            self.add_linear(i, val)
        else:
            key = (i, j) if i < j else (j, i)
            self.Q[key] = self.Q.get(key, 0.0) + val

    # ---------------------------------

    def add_one_hot(self):
        """
        One-Hot Constraint: P*(sum_i x_i - 1)^2
        = P*(sum x_i^2 + 2*sum_{i<j} x_i*x_j - 2*sum x_i + 1)
        Binär: x_i^2 = x_i  →  Diagonal: -P, Kreuzterme: +2P
        """

        for var in self.domains:

            idx_list = list(
                self.var_map[var].values()
            )

            txt = "+".join(
                f"{var}_{v}"
                for v in self.domains[var]
            )

            self.energy_terms.append(
                f"{self.P}*({txt}-1)^2"
            )

            for i in idx_list:
                self.add_linear(i, -self.P)

            for i, j in combinations(idx_list, 2):
                self.add_quad(i, j, 2 * self.P)

    # ---------------------------------

    def translate_literal(self, expr):
        """
        Übersetzt ein Literal in (index_liste, text).

        x==v  →  direkt ein einzelner Binärindex
        x<=v, x<v, x>=v, x>v  →  Summe mehrerer Vars → Hilfsvariable h
                  via P*(h-sum)^2. Dank One-Hot gilt sum ∈ {0,1},
                  also h ∈ {0,1} korrekt.
        """

        # Reihenfolge wichtig: <=/>= muessen vor </> geprueft werden,
        # sonst matcht "<=" faelschlich bereits als "<".
        m = re.match(
            r"([a-zA-Z]\w*)(==|<=|>=|<|>)(\d+)",
            expr.strip()
        )

        if not m:
            raise Exception(f"Unbekanntes Literal: {expr}")

        var   = m.group(1)
        op    = m.group(2)
        value = int(m.group(3))

        if var not in self.domains:
            raise Exception(
                f"Variable '{var}' wurde nicht deklariert. "
                f"Bitte bei den Wertebereichen mit angeben "
                f"(z.B. '{var}=1-5')."
            )

        if op == "==":
            if value not in self.var_map[var]:
                raise Exception(
                    f"Wert {value} liegt außerhalb des Wertebereichs "
                    f"von '{var}' ({min(self.domains[var])}-"
                    f"{max(self.domains[var])})."
                )
            idx = self.var_map[var][value]
            return [idx], f"{var}_{value}"

        # Vergleichsoperator generisch anwenden
        predicates = {
            "<=": lambda v: v <= value,
            ">=": lambda v: v >= value,
            "<":  lambda v: v <  value,
            ">":  lambda v: v >  value,
        }
        match_fn = predicates[op]

        sub_idxs = [
            self.var_map[var][v]
            for v in self.domains[var]
            if match_fn(v)
        ]

        sub_txt = "+".join(
            f"{var}_{v}"
            for v in self.domains[var]
            if match_fn(v)
        )

        # Kein Domänenwert erfüllt das Literal (z.B. x>10 bei Domäne 1-5)
        # → h muss immer 0 sein, kein Summenterm nötig.
        if len(sub_idxs) == 0:
            h_idx, h_name = self.create_aux("h")
            self.add_linear(h_idx, self.P)
            self.energy_terms.append(f"{self.P}*{h_name}")
            return [h_idx], h_name

        # Nur ein Wert passt → keine Hilfsvariable nötig
        if len(sub_idxs) == 1:
            return sub_idxs, sub_txt

        # Hilfsvariable h
        h_idx, h_name = self.create_aux("h")

        self.energy_terms.append(
            f"{self.P}*({h_name}-({sub_txt}))^2"
        )

        P = self.P

        # P*(h - sum)^2 ausmultipliziert:
        # h^2 = h (binär) → +P auf Diagonale h
        self.add_linear(h_idx, P)

        # -2*h*x_i
        for i in sub_idxs:
            self.add_quad(h_idx, i, -2 * P)

        # sum^2 = sum_i x_i + 2*sum_{i<j} x_i*x_j  (binär)
        for i in sub_idxs:
            self.add_linear(i, P)

        for i, j in combinations(sub_idxs, 2):
            self.add_quad(i, j, 2 * P)

        return [h_idx], h_name

    # ---------------------------------

    def add_or_two(
        self,
        left_idxs,
        right_idxs,
        left_txt,
        right_txt
    ) :
        """
        Kodiert z = OR(L, R) mit Hilfsvariable z und Slack u.
        L, R sind einzelne Binärvariablen ∈ {0,1}.

        Penalty 1: P*(z - L - R + u)^2   →  erzwingt z = OR(L,R)
        Penalty 2: P*(L*R - 2u*L - 2u*R + 3u)  →  erzwingt u = AND(L,R)
        """

        z_idx, z_name = self.create_aux("z")
        u_idx, u_name = self.create_aux("u")

        P = self.P

        self.energy_terms.append(
            f"{P}*({z_name}-({left_txt})-({right_txt})+{u_name})^2"
        )
        self.energy_terms.append(
            f"{P}*(({left_txt})*({right_txt})"
            f"-2*{u_name}*({left_txt})"
            f"-2*{u_name}*({right_txt})"
            f"+3*{u_name})"
        )

        # --- Penalty 1: P*(z - L - R + u)^2 ---
        self.add_linear(z_idx, P)       # z^2 = z
        self.add_linear(u_idx, P)       # u^2 = u

        for i in left_idxs:
            self.add_linear(i, P)       # L^2 = L
        for i in right_idxs:
            self.add_linear(i, P)       # R^2 = R

        for i in left_idxs:
            self.add_quad(z_idx, i, -2 * P)   # -2*z*L
        for i in right_idxs:
            self.add_quad(z_idx, i, -2 * P)   # -2*z*R

        self.add_quad(z_idx, u_idx, 2 * P)    # +2*z*u

        for i in left_idxs:
            for j in right_idxs:
                self.add_quad(i, j, 2 * P)    # +2*L*R

        for i in left_idxs:
            self.add_quad(u_idx, i, -2 * P)   # -2*u*L
        for i in right_idxs:
            self.add_quad(u_idx, i, -2 * P)   # -2*u*R

        # --- Penalty 2: P*(L*R - 2u*L - 2u*R + 3u) ---
        for i in left_idxs:
            for j in right_idxs:
                self.add_quad(i, j, P)         # +L*R

        for i in left_idxs:
            self.add_quad(u_idx, i, -2 * P)   # -2*u*L
        for i in right_idxs:
            self.add_quad(u_idx, i, -2 * P)   # -2*u*R

        self.add_linear(u_idx, 3 * P)          # +3*u

        return [z_idx], z_name

    # ---------------------------------

    def parse_formula(self, formula):
        """
        Parst KNF: (L1 OR L2 OR ...) AND (L3 OR L4) AND ...
        Beliebig viele Literale pro Klausel durch iterative OR-Reduktion.
        """

        clauses = re.findall(
            r"\((.*?)\)",
            formula
        )

        for clause in clauses:

            # \bOR\b statt split("OR"), damit Variablennamen,
            # die "OR" als Teilstring enthalten (z.B. "door"), nicht
            # faelschlich zerteilt werden.
            literals = [
                l.strip()
                for l in re.split(r"\bOR\b", clause)
            ]

            if len(literals) == 1:
                # Einzelliteral: P*(1-L)  →  -P auf Diagonale
                idxs, txt = self.translate_literal(literals[0])
                self.energy_terms.append(
                    f"{self.P}*(1-({txt}))"
                )
                for i in idxs:
                    self.add_linear(i, -self.P)

            else:
                # Erstes OR
                left_idxs,  left_txt  = self.translate_literal(literals[0])
                right_idxs, right_txt = self.translate_literal(literals[1])

                current_idxs, current_txt = self.add_or_two(
                    left_idxs, right_idxs,
                    left_txt,  right_txt
                )

                # Weitere Literale iterativ verknüpfen
                for lit in literals[2:]:
                    next_idxs, next_txt = self.translate_literal(lit)
                    current_idxs, current_txt = self.add_or_two(
                        current_idxs, next_idxs,
                        current_txt,  next_txt
                    )

                # Finale Klausel-Penalty: P*(1-z)  →  -P auf z-Diagonale
                self.energy_terms.append(
                    f"{self.P}*(1-{current_txt})"
                )
                for i in current_idxs:
                    self.add_linear(i, -self.P)

    # ---------------------------------

    def solve(self, Q, num_reads=1000):
        """
        Löst das QUBO mit Simulated Annealing (dwave-samplers).
        Q ist bereits ein Dictionary {(i, j): wert} — kein Umwandeln
        aus einer dichten Matrix mehr nötig.
        """

        print("\n" + "=" * 50)
        print("SOLVER  —  Simulated Annealing")
        print("=" * 50)

        real = sum(
            len(v) for v in self.domains.values()
        )

        print(
            f"Variablen : {self.n} "
            f"({real} echte + {self.n - real} Hilfsvariablen)"
        )
        print(f"num_reads : {num_reads}")

        Q_dict = {k: v for k, v in Q.items() if v != 0.0}

        sampler  = SimulatedAnnealingSampler()
        response = sampler.sample_qubo(
            Q_dict,
            num_reads=num_reads,
            num_sweeps=10000
        )

        best = response.first.sample

        print(f"\nBeste Energie : {response.first.energy:.4f}")

        # Alle Binärvariablen ausgeben
        print("\n--- Belegung aller Binärvariablen ---")
        for i, name in enumerate(self.variables):
            val    = best.get(i, 0)
            marker = "  ←" if val == 1 else ""
            print(f"  {i:3d}  {name:<14s} = {val}{marker}")

        # One-Hot dekodieren
        print("\n--- Lösung (One-Hot dekodiert) ---")
        violations  = []
        assignment  = {}

        for var in self.domains:
            chosen = [
                v for v in self.domains[var]
                if best.get(self.var_map[var][v], 0) == 1
            ]
            if len(chosen) == 1:
                print(f"  {var} = {chosen[0]}")
                assignment[var] = chosen[0]
            else:
                violations.append(var)
                print(f"  {var} = ??? (One-Hot verletzt, gewählt: {chosen})")

        print()
        if not violations:
            print("✓  Alle One-Hot Constraints erfüllt.")
        else:
            print(f"⚠  One-Hot verletzt für: {', '.join(violations)}")
            print("   → Penalty erhöhen oder num_reads vergrößern.")

        # Formel mit der gefundenen Belegung auswerten
        print("\n--- SAT Formel (substituiert) ---")
        substituted, overall_result = self.evaluate_formula(assignment)
        print(f"  {substituted}")
        print(f"  = {'TRUE' if overall_result else 'FALSE'}")

        return response

    # ---------------------------------

    def print_labeled_matrix(self, M, decimals=1, output_file=None):
        """
        Druckt die QUBO-Matrix mit Zeilen-/Spaltenbeschriftung
        (Variablennamen statt reiner Indizes).

        Hinweis: bei vielen Variablen (viele Hilfsvariablen z/u/h)
        wird die Ausgabe entsprechend breit -- im Terminal kann das
        zu automatischem Zeilenumbruch führen. Über output_file kann
        die vollständige, unumgebrochene Ausgabe zusätzlich in eine
        Datei geschrieben werden (z.B. output.txt), die sich ohne
        Zeilenumbruch öffnen lässt.
        """

        names = self.variables

        # Spaltenbreite an den längsten Namen anpassen
        width = max(max(len(n) for n in names), 6) + 2

        lines = []

        header = " " * width + "".join(
            f"{n:>{width}}" for n in names
        )
        lines.append(header)

        for i, row in enumerate(M):
            line = f"{names[i]:<{width}}" + "".join(
                f"{val:>{width}.{decimals}f}" for val in row
            )
            lines.append(line)

        text = "\n".join(lines)

        print(text)

        if output_file:
            with open(output_file, "w") as f:
                f.write(text + "\n")

    # ---------------------------------

    def evaluate_formula(self, assignment):
        """
        Setzt die gefundene Variablenbelegung (var -> Wert) in die
        ursprüngliche SAT-Formel ein und gibt sowohl den mit
        TRUE/FALSE substituierten Ausdruck als auch das
        Gesamtergebnis (True/False) zurück.

        Wichtig: Die Auswertung erfolgt anhand der echten,
        one-hot-dekodierten Variablenwerte (assignment), nicht
        anhand der QUBO-Hilfsvariablen h/z/u — dadurch ist das
        Ergebnis auch dann semantisch korrekt, wenn Hilfsvariablen
        im Sample zufällig "falsch" gelandet sind.
        """

        predicates = {
            "==": lambda a, b: a == b,
            "<=": lambda a, b: a <= b,
            ">=": lambda a, b: a >= b,
            "<":  lambda a, b: a <  b,
            ">":  lambda a, b: a >  b,
        }

        clauses = re.findall(r"\((.*?)\)", self.formula)

        clause_strs   = []
        overall_result = True

        for clause in clauses:

            literals = [
                l.strip()
                for l in re.split(r"\bOR\b", clause)
            ]

            lit_labels = []
            clause_result = False

            for lit in literals:

                m = re.match(
                    r"([a-zA-Z]\w*)(==|<=|>=|<|>)(\d+)",
                    lit.strip()
                )
                var, op, value = m.group(1), m.group(2), int(m.group(3))

                var_value = assignment.get(var)

                # Falls die Variable (z.B. wegen One-Hot-Verletzung)
                # keinen eindeutigen Wert hat, gilt das Literal als
                # nicht auswertbar / falsch.
                if var_value is None:
                    lit_result = False
                else:
                    lit_result = predicates[op](var_value, value)

                lit_labels.append("TRUE" if lit_result else "FALSE")
                clause_result = clause_result or lit_result

            clause_strs.append(
                "(" + " OR ".join(lit_labels) + ")"
            )
            overall_result = overall_result and clause_result

        substituted = " AND ".join(clause_strs)

        return substituted, overall_result

    # ---------------------------------

    def to_dense_matrix(self):
        """
        Erzeugt aus dem Q-Dictionary eine dichte Matrix der Größe
        self.n x self.n (exakt passend, nicht mehr hartkodiert).
        Nur zur Anzeige/Inspektion gedacht — der Solver arbeitet
        weiterhin direkt mit dem Dictionary.
        """

        M = np.zeros((self.n, self.n))

        for (i, j), val in self.Q.items():
            M[i][j] = val

        return M

    # ---------------------------------

    def build(self):

        self.add_one_hot()

        formula = self.read_formula()
        self.formula = formula

        self.parse_formula(formula)

        return self.Q


# =======================

if __name__ == "__main__":

    compiler = SATtoQUBO()

    compiler.read_domains()

    Q = compiler.build()

    print("\n" + "=" * 50)
    print("ENERGIEFUNKTION")
    print("=" * 50)
    print("H(x) =")
    for t in compiler.energy_terms:
        print("  +", t)

    print("\n" + "=" * 50)
    print("QUBO MATRIX")
    print("=" * 50)
    print(f"Anzahl Variablen : {compiler.n}")
    M = compiler.to_dense_matrix()
    compiler.print_labeled_matrix(M, output_file="output.txt")

    print("\n" + "=" * 50)
    print("VARIABLEN")
    print("=" * 50)
    for i, v in enumerate(compiler.variables):
        print(f"  {i:3d} : {v}")

    compiler.solve(Q, num_reads=1000)