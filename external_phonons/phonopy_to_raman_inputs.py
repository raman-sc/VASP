import argparse
import json
import math
from pathlib import Path

import yaml
from ase.io import read
from ase.data import atomic_masses, atomic_numbers

THZ_TO_CM = 33.35640951981521

parser = argparse.ArgumentParser(
    description="Convert Phonopy qpoints.yaml to raman-sc/VASP freq.dat and modes_sqrt_amu.dat"
)
parser.add_argument("--qpoints", required=True, help="Path to qpoints.yaml")
parser.add_argument("--poscar", required=True, help="Path to POSCAR.phon or matching POSCAR")
parser.add_argument("--outdir", required=True, help="Output directory")
args = parser.parse_args()

qpoints_path = Path(args.qpoints)
poscar_path = Path(args.poscar)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

atoms = read(poscar_path, format="vasp")
symbols = atoms.get_chemical_symbols()
masses = [atomic_masses[atomic_numbers[s]] for s in symbols]
nat = len(symbols)

with qpoints_path.open() as f:
    data = yaml.safe_load(f)

qblock = data["phonon"][0]
bands = qblock["band"]

if len(bands) != 3 * nat:
    raise SystemExit(f"Mode number mismatch: {len(bands)} modes but 3N = {3 * nat}")

modes = []
for i, band in enumerate(bands, start=1):
    freq_thz = float(band["frequency"])
    freq_cm = freq_thz * THZ_TO_CM

    eig = []
    for atom_vec in band["eigenvector"]:
        v = []
        for comp in atom_vec:
            if isinstance(comp, list):
                v.append(float(comp[0]))
            else:
                v.append(float(comp))
        eig.append(v)

    if len(eig) != nat:
        raise SystemExit(f"Atom number mismatch in mode {i}: {len(eig)} vs {nat}")

    modes.append({
        "original_mode": i,
        "frequency_THz": freq_thz,
        "frequency_cm-1": freq_cm,
        "eigenvector": eig,
    })

acoustic = [m for m in modes if abs(m["frequency_cm-1"]) < 1.0]
optical = [m for m in modes if abs(m["frequency_cm-1"]) >= 1.0]

ordered = sorted(optical, key=lambda m: m["frequency_cm-1"], reverse=True)
ordered += sorted(acoustic, key=lambda m: m["frequency_cm-1"], reverse=True)

with (outdir / "freq.dat").open("w") as f:
    for m in ordered:
        f.write(f'{m["frequency_cm-1"]:.10f}\n')

with (outdir / "modes_sqrt_amu.dat").open("w") as f:
    for m in ordered:
        for atom_i, vec in enumerate(m["eigenvector"]):
            sqrt_mass = math.sqrt(masses[atom_i])
            f.write(
                f"{vec[0]/sqrt_mass: .12e} "
                f"{vec[1]/sqrt_mass: .12e} "
                f"{vec[2]/sqrt_mass: .12e}\n"
            )
        f.write("\n")

with (outdir / "phonopy_to_raman_mode_order.csv").open("w") as f:
    f.write("new_mode,original_phonopy_mode,frequency_cm-1,frequency_THz\n")
    for new_i, m in enumerate(ordered, start=1):
        f.write(
            f'{new_i},{m["original_mode"]},'
            f'{m["frequency_cm-1"]:.10f},{m["frequency_THz"]:.10f}\n'
        )

with (outdir / "conversion_log.txt").open("w") as f:
    f.write("Phonopy to raman-sc/VASP conversion log\n")
    f.write(f"qpoints file: {qpoints_path}\n")
    f.write(f"POSCAR file: {poscar_path}\n")
    f.write(f"Number of atoms: {nat}\n")
    f.write(f"Number of modes: {len(ordered)}\n")
    f.write("Mass scaling: eigenvectors divided by sqrt(atomic mass in amu)\n")
    f.write("Atom symbols and masses:\n")
    for s, m in zip(symbols, masses):
        f.write(f"  {s}: {m}\n")

print("Wrote:")
print(outdir / "freq.dat")
print(outdir / "modes_sqrt_amu.dat")
print(outdir / "phonopy_to_raman_mode_order.csv")
print(outdir / "conversion_log.txt")
print()
print("First 5 converted modes:")
for i, m in enumerate(ordered[:5], start=1):
    print(i, m["original_mode"], f'{m["frequency_cm-1"]:.6f}', "cm^-1")
