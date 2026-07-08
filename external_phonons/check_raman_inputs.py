import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Check external phonon input files for raman-sc/VASP."
)
parser.add_argument("--poscar", required=True)
parser.add_argument("--freq", required=True)
parser.add_argument("--modes", required=True)
args = parser.parse_args()

poscar = Path(args.poscar)
freq = Path(args.freq)
modes = Path(args.modes)

def read_poscar_natoms(path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    line5 = lines[5].split()

    def all_ints(items):
        try:
            [int(x) for x in items]
            return True
        except ValueError:
            return False

    if all_ints(line5):
        counts = [int(x) for x in line5]
    else:
        counts = [int(x) for x in lines[6].split()]

    return sum(counts)

errors = []

if not poscar.exists():
    errors.append(f"Missing POSCAR file: {poscar}")
if not freq.exists():
    errors.append(f"Missing freq.dat file: {freq}")
if not modes.exists():
    errors.append(f"Missing modes_sqrt_amu.dat file: {modes}")

if errors:
    for e in errors:
        print("ERROR:", e)
    raise SystemExit(1)

natoms = read_poscar_natoms(poscar)
expected_modes = 3 * natoms

freq_lines = [line for line in freq.read_text().splitlines() if line.strip()]

blocks = []
current = []
for line in modes.read_text().splitlines():
    if not line.strip():
        if current:
            blocks.append(current)
            current = []
    else:
        parts = line.split()
        if len(parts) != 3:
            errors.append(f"Bad mode line: {line}")
        current.append(parts)

if current:
    blocks.append(current)

if len(freq_lines) != expected_modes:
    errors.append(f"freq.dat has {len(freq_lines)} modes, expected {expected_modes}")

if len(blocks) != expected_modes:
    errors.append(f"modes_sqrt_amu.dat has {len(blocks)} mode blocks, expected {expected_modes}")

for i, block in enumerate(blocks, start=1):
    if len(block) != natoms:
        errors.append(f"Mode {i} has {len(block)} atom rows, expected {natoms}")

print(f"Number of atoms: {natoms}")
print(f"Expected modes: {expected_modes}")
print(f"freq.dat modes: {len(freq_lines)}")
print(f"modes_sqrt_amu.dat blocks: {len(blocks)}")

if errors:
    print("Status: FAIL")
    for e in errors:
        print("ERROR:", e)
    raise SystemExit(1)

print("Status: PASS")
