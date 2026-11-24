#!/usr/bin/env python3
"""General .reaction file writer for Ancalagon (any beam/target, any
branching level scheme, optional BKIN/RWID/RTUN), for the DRAGON NN
emulator campaigns. generate_reactions.py stays as the k39-specific,
100%-linear-chain generator it always was; this module is the general
version.

Card format: see include/ReactionConfig.hh.

Terminus shift (the "RECL trick" from generate_reactions.py, generalised):
Ancalagon always ends a cascade at level 0, the recoil's RECL mass. Writing
RECL = recoil_mass_excess + E_terminus makes "level 0" a chosen excited
state (e.g. 40Ca 0+_2, whose E0 decay Ancalagon cannot emit) or an
arbitrary fictional energy; LEVL energies are then excitation above that
terminus. E_terminus = 0 is the true ground state.

Usage as a library:

    from reaction_writer import NUCLEI, Level, Branch, ReactionSpec
    spec = ReactionSpec(beam=NUCLEI["39K"], target=NUCLEI["1H"],
                        recoil=NUCLEI["40Ca"], charge_state=8, eres_mev=0.606,
                        levels={1: Level(3.73669)}, branches=[Branch(-1, 80, 0),
                        Branch(-1, 20, 1), Branch(1, 100, 0)])
    spec.write("out.reaction")

or from the command line with a JSON spec (same field names):

    python3 reaction_writer.py spec.json out.reaction
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field

AMU_MEV = 931.49432  # same constant as src/ReactionKinematics.cc


@dataclass(frozen=True)
class Nucleus:
    Z: int
    A: int
    mass_excess_mev: float
    label: str

    @property
    def mass_mev(self) -> float:
        return self.A * AMU_MEV + self.mass_excess_mev


# AME2020 mass excesses (MeV) for the validation reactions in the NN
# emulator plan. Check against the AME table before relying on a new entry.
# 1H, 39K and 40Ca reproduce the values in the existing k39 files.
NUCLEI = {
    "1H": Nucleus(1, 1, 7.28897, "1H"),
    "4He": Nucleus(2, 4, 2.42492, "4He"),
    "15O": Nucleus(8, 15, 2.85561, "15O"),
    "19F": Nucleus(9, 19, -1.48744, "19F"),
    "19Ne": Nucleus(10, 19, 1.75205, "19Ne"),
    "20Ne": Nucleus(10, 20, -7.04193, "20Ne"),
    "22Ne": Nucleus(10, 22, -8.02472, "22Ne"),
    "23Na": Nucleus(11, 23, -9.52985, "23Na"),
    "26Mg": Nucleus(12, 26, -16.21455, "26Mg"),
    "38K": Nucleus(19, 38, -28.80075, "38K"),
    "39Ca": Nucleus(20, 39, -27.28271, "39Ca"),
    "39K": Nucleus(19, 39, -33.8072, "39K"),
    "40Ca": Nucleus(20, 40, -34.8460, "40Ca"),
}


@dataclass
class Level:
    excitation_mev: float          # above the terminus (see module docstring)
    lifetime_s: float = 1.0e-12    # parsed but unused by Ancalagon kinematics


@dataclass
class Branch:
    from_level: int   # -1 = resonance
    percent: float
    to_level: int     # 0 = terminus


@dataclass
class ReactionSpec:
    beam: Nucleus
    target: Nucleus
    recoil: Nucleus
    charge_state: int
    eres_mev: float
    levels: dict[int, Level] = field(default_factory=dict)
    branches: list[Branch] = field(default_factory=list)
    terminus_mev: float = 0.0       # excitation of the state used as level 0
    rwid_mev: float | None = None   # RWID: CM Breit-Wigner width
    bkin_mev: float | None = None   # BKIN: beam lab KE entering the gas
    rtun: float | None = None       # RTUN: fixed separator retune scale
    comments: list[str] = field(default_factory=list)

    # --- kinematics helpers -------------------------------------------
    @property
    def q_value_mev(self) -> float:
        return self.beam.mass_excess_mev + self.target.mass_excess_mev - self.recoil.mass_excess_mev

    @property
    def ex_resonance_mev(self) -> float:
        """Resonance excitation above the recoil's TRUE ground state."""
        return self.q_value_mev + self.eres_mev

    @property
    def e_budget_mev(self) -> float:
        """Energy shared by the cascade gammas (resonance -> terminus)."""
        return self.ex_resonance_mev - self.terminus_mev

    @property
    def recl_mass_excess_mev(self) -> float:
        return self.recoil.mass_excess_mev + self.terminus_mev

    def beam_ke_for_ecm(self, ecm_mev: float) -> float:
        """Lab beam KE (MeV) for a CM energy above threshold, target at rest
        (same formula as ReactionKinematics::BeamKineticEnergyMeVForEcm)."""
        mb, mt = self.beam.mass_mev, self.target.mass_mev
        W = mb + mt + ecm_mev
        return (W * W - mb * mb - mt * mt) / (2.0 * mt) - mb

    @property
    def beam_ke_resonance_mev(self) -> float:
        return self.beam_ke_for_ecm(self.eres_mev)

    def compound_beta(self, ecm_mev: float | None = None) -> float:
        ke = self.beam_ke_for_ecm(self.eres_mev if ecm_mev is None else ecm_mev)
        e_tot = ke + self.beam.mass_mev + self.target.mass_mev
        p = math.sqrt((ke + self.beam.mass_mev) ** 2 - self.beam.mass_mev ** 2)
        return p / e_tot

    # --- validation / cascade enumeration --------------------------------
    def paths(self) -> list[tuple[float, list[int]]]:
        """Every resonance -> terminus path as (probability, [levels]),
        raising on cycles, dead ends or undefined levels."""
        by_from: dict[int, list[Branch]] = {}
        for b in self.branches:
            if b.percent <= 0:
                raise ValueError(f"branch {b} has non-positive percent")
            if b.from_level != -1 and b.from_level not in self.levels:
                raise ValueError(f"branch from undefined level {b.from_level}")
            if b.to_level != 0 and b.to_level not in self.levels:
                raise ValueError(f"branch to undefined level {b.to_level}")
            by_from.setdefault(b.from_level, []).append(b)

        out: list[tuple[float, list[int]]] = []

        def walk(level: int, prob: float, seen: list[int]) -> None:
            if level == 0:
                out.append((prob, seen + [0]))
                return
            if level in seen:
                raise ValueError(f"cascade cycle through level {level}")
            opts = by_from.get(level)
            if not opts:
                raise ValueError(f"level {level} has no outgoing BRAT -- dead end")
            tot = sum(b.percent for b in opts)
            for b in opts:
                walk(b.to_level, prob * b.percent / tot, seen + [level])

        walk(-1, 1.0, [])
        return out

    def level_energy(self, level: int) -> float:
        if level == -1:
            return self.e_budget_mev
        if level == 0:
            return 0.0
        return self.levels[level].excitation_mev

    def gamma_energies(self, path: list[int]) -> list[float]:
        """Transition energies (MeV, before recoil corrections) along a path."""
        return [self.level_energy(a) - self.level_energy(b) for a, b in zip(path[:-1], path[1:])]

    def validate(self) -> None:
        if self.e_budget_mev <= 0:
            raise ValueError("terminus lies above the resonance")
        for n, lv in self.levels.items():
            if n <= 0:
                raise ValueError("LEVL numbers must be positive")
            if not (0.0 < lv.excitation_mev < self.e_budget_mev):
                raise ValueError(f"level {n} at {lv.excitation_mev} MeV is outside (0, e_budget)")
        for b in self.branches:
            if self.level_energy(b.from_level) <= self.level_energy(b.to_level):
                raise ValueError(f"branch {b} does not go down in energy")
        self.paths()

    # --- output ----------------------------------------------------------
    def text(self) -> str:
        self.validate()
        L = []
        for c in self.comments:
            L.append(f"COMM {c}")
        L.append(f"COMM {self.beam.label}({self.target.label},g){self.recoil.label}: "
                 f"Q={self.q_value_mev:.5f} MeV, Ex(res)={self.ex_resonance_mev:.5f} MeV, "
                 f"terminus={self.terminus_mev:.5f} MeV, e_budget={self.e_budget_mev:.5f} MeV")
        L.append(f"COMM Beam KE at resonance {self.beam_ke_resonance_mev:.5f} MeV, "
                 f"compound beta {self.compound_beta():.6f}. Written by reaction_writer.py.")
        for p, path in self.paths():
            g = ", ".join(f"{e:.4f}" for e in self.gamma_energies(path))
            L.append(f"COMM path {path} prob {p:.5f}: gammas [{g}] MeV")
        L.append("")
        if self.rtun is not None:
            L.append(f"RTUN {self.rtun:.6f}")
        L.append(f"BEAM {self.beam.Z:<3d} {self.beam.A:<4d} {self.beam.mass_excess_mev:.5f}")
        L.append(f"TARG {self.target.Z:<3d} {self.target.A:<4d} {self.target.mass_excess_mev:.5f}")
        L.append(f"RECL {self.recoil.Z:<3d} {self.recoil.A:<4d} {self.recl_mass_excess_mev:.5f} "
                 f"{self.charge_state}")
        L.append(f"ERES {self.eres_mev:.6f}")
        if self.rwid_mev is not None:
            L.append(f"RWID {self.rwid_mev:.6e}")
        if self.bkin_mev is not None:
            L.append(f"BKIN {self.bkin_mev:.5f}")
        L.append("")
        for n in sorted(self.levels):
            lv = self.levels[n]
            L.append(f"LEVL {n:<3d} {lv.excitation_mev:.5f} {lv.lifetime_s:.3E}")
        L.append("")
        for b in self.branches:
            L.append(f"BRAT {b.from_level:<3d} {b.percent:.6f} {b.to_level}")
        L.append("SENT")
        return "\n".join(L) + "\n"

    def write(self, path: str) -> None:
        with open(path, "w") as f:
            f.write(self.text())

    # --- JSON ------------------------------------------------------------
    @classmethod
    def from_dict(cls, d: dict) -> "ReactionSpec":
        def nuc(x):
            return NUCLEI[x] if isinstance(x, str) else Nucleus(**x)
        return cls(
            beam=nuc(d["beam"]), target=nuc(d["target"]), recoil=nuc(d["recoil"]),
            charge_state=int(d["charge_state"]), eres_mev=float(d["eres_mev"]),
            levels={int(k): Level(**v) if isinstance(v, dict) else Level(float(v))
                    for k, v in d.get("levels", {}).items()},
            branches=[Branch(*b) if isinstance(b, (list, tuple)) else Branch(**b)
                      for b in d.get("branches", [])],
            terminus_mev=float(d.get("terminus_mev", 0.0)),
            rwid_mev=d.get("rwid_mev"), bkin_mev=d.get("bkin_mev"), rtun=d.get("rtun"),
            comments=list(d.get("comments", [])),
        )


def linear_chain(spec_base: dict, gammas_mev: list[float]) -> ReactionSpec:
    """A 100% linear cascade emitting `gammas_mev` (top first). The gammas
    must sum to the e_budget implied by spec_base (beam/target/recoil/eres/
    terminus); the last one is adjusted for rounding and the rest checked."""
    spec = ReactionSpec.from_dict({**spec_base, "levels": {}, "branches": []})
    budget = spec.e_budget_mev
    if abs(sum(gammas_mev) - budget) > 1e-3:
        raise ValueError(f"gammas sum {sum(gammas_mev):.5f} != e_budget {budget:.5f}")
    n = len(gammas_mev)
    levels, e = {}, budget
    for i, g in enumerate(gammas_mev[:-1]):
        e -= g
        levels[n - 1 - i] = Level(e)
    spec.levels = levels
    if n == 1:
        spec.branches = [Branch(-1, 100.0, 0)]
    else:
        spec.branches = [Branch(-1, 100.0, n - 1)] + [Branch(k, 100.0, k - 1) for k in range(n - 1, 0, -1)]
    return spec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("spec_json")
    ap.add_argument("out")
    a = ap.parse_args()
    with open(a.spec_json) as f:
        ReactionSpec.from_dict(json.load(f)).write(a.out)


if __name__ == "__main__":
    main()
