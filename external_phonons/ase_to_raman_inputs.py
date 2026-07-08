import json
from pathlib import Path

INPUT = Path("ase_gamma_modes.json")
OUTDIR = Path("ase_raman_inputs")
OUTDIR.mkdir(exist_ok=True)

with INPUT.open() as f:
    modes = json.load(f)

# Remove imaginary/acoustic modes from the front only for ordering.
acoustic = [m for m in modes if abs(m["frequency_cm-1"]) < 1.0]
optical = [m for m in modes if abs(m["frequency_cm-1"]) >= 1.0]

# Match raman-sc ordering: high-frequency optical modes first, acoustic modes last.
ordered = sorted(optical, key=lambda m: m["frequency_cm-1"], reverse=True)
ordered += sorted(acoustic, key=lambda m: m["frequency_cm-1"], reverse=True)

with (OUTDIR / "freq.dat").open("w") as f:
    for m in ordered:
        f.write(f'{m["frequency_cm-1"]:.10f}\n')

with (OUTDIR / "modes_sqrt_amu.dat").open("w") as f:
    for m in ordered:
        for v in m["eigenvector"]:
            f.write(f"{v[0]: .12e} {v[1]: .12e} {v[2]: .12e}\n")
        f.write("\n")

with (OUTDIR / "ase_to_raman_mode_order.csv").open("w") as f:
    f.write("new_mode,original_ase_mode,frequency_cm-1\n")
    for i, m in enumerate(ordered, start=1):
        f.write(f'{i},{m["mode"]},{m["frequency_cm-1"]:.10f}\n')

print("Wrote ase_raman_inputs/freq.dat")
print("Wrote ase_raman_inputs/modes_sqrt_amu.dat")
print("Wrote ase_raman_inputs/ase_to_raman_mode_order.csv")
print("First 6 modes:")
for i, m in enumerate(ordered[:6], start=1):
    print(i, m["mode"], f'{m["frequency_cm-1"]:.6f}')
