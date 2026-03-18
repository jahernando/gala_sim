#!/bin/bash
# Run the Geant4 photon simulation with correct environment
# Usage: ./run_geant4.sh [input.json] [output.root]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="${1:-$SCRIPT_DIR/el_electrons.json}"
OUTPUT="${2:-$SCRIPT_DIR/build/output.root}"

# Geant4 data paths (conda-forge layout)
DATADIR="${CONDA_PREFIX}/share/Geant4/data"
export G4NEUTRONHPDATA="$DATADIR/NDL4.6"
export G4LEDATA="$DATADIR/EMLOW8.0"
export G4LEVELGAMMADATA="$DATADIR/PhotonEvaporation5.7"
export G4RADIOACTIVEDATA="$DATADIR/RadioactiveDecay5.6"
export G4PARTICLEXSDATA="$DATADIR/PARTICLEXS4.0"
export G4PIIDATA="$DATADIR/PII1.3"
export G4REALSURFACEDATA="$DATADIR/RealSurface2.2"
export G4SAIDXSDATA="$DATADIR/SAIDDATA2.0"
export G4ABLADATA="$DATADIR/ABLA3.1"
export G4INCLDATA="$DATADIR/INCL1.0"
export G4ENSDFSTATEDATA="$DATADIR/ENSDFSTATE2.3"

echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo "Running ELHolePhotonSim..."

"$SCRIPT_DIR/build/ELHolePhotonSim" "$INPUT" "$OUTPUT"
