// ============================================================================
//  ELHolePhotonSim.cc  — Geant4 VUV photon propagation in NEXT honeycomb hole
// ============================================================================
//
//  PURPOSE:
//    Reads el_electrons.json (excitation points from Python simulation),
//    generates VUV photons at each excitation point, propagates them through
//    the cylindrical hole geometry, and records:
//      - arrival position at SiPM (hole bottom)
//      - number of wall reflections
//      - whether photon was absorbed or escaped
//
//  GEOMETRY (single hole, centred at origin):
//
//    z = HOLE_HEIGHT (0.5 cm)  ─── hole top (open to drift region)
//    │                         │
//    │   Xe gas (EL region)    │   ← excitation points here
//    │   R = 2.5 mm            │
//    │   PTFE cylindrical wall │
//    │                         │
//    z = 0                     ─── SiPM (hole bottom)
//
//  PHYSICS:
//    - VUV photon λ = 172 nm (Xe scintillation)
//    - PTFE wall reflectivity: 95% (diffuse Lambertian)
//    - SiPM: perfect absorber (detector)
//    - Photons launched isotropically from excitation points
//
//  BUILD:
//    mkdir build && cd build
//    cmake .. -DGeant4_DIR=/path/to/geant4/lib/Geant4-11.X.X
//    make -j4
//    ./ELHolePhotonSim ../el_electrons.json output.root
//
//  REQUIRES: Geant4 >= 10.7, nlohmann/json (auto-fetched by CMake)
// ============================================================================

#include "G4RunManager.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4ParticleGun.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────
//  Geometry constants (match Python simulation)
// ─────────────────────────────────────────────────────────────────
static const G4double HOLE_RADIUS   = 0.25  * cm;   // 2.5 mm
static const G4double HOLE_HEIGHT   = 0.50  * cm;   // 5 mm (EL gap)
static const G4double WALL_THICK    = 0.25  * cm;   // wall around hole
static const G4double WORLD_HALF    = 2.0   * cm;
static const G4double PTFE_REFLECT  = 0.95;          // VUV reflectivity
static const G4double XE_WL_NM     = 172.0;          // Xe scintillation [nm]

// ─────────────────────────────────────────────────────────────────
//  Electron data loaded from JSON
// ─────────────────────────────────────────────────────────────────
struct ExcitationPoint {
    G4double x, y, z;   // in G4 internal units (mm)
    G4int    n_photons;
};

struct ElectronRecord {
    G4int    id;
    std::vector<ExcitationPoint> excitations;
    G4int    total_photons;
};

std::vector<ElectronRecord> gElectrons;

void LoadElectronData(const std::string& fname) {
    std::ifstream f(fname);
    if (!f) throw std::runtime_error("Cannot open: " + fname);
    json j; f >> j;

    for (auto& e : j["electrons"]) {
        // Only load electrons that entered the hole and produced photons
        if (!e["entered_hole"].get<bool>()) continue;
        if (e["excitation_points"].empty()) continue;

        ElectronRecord rec;
        rec.id = e["electron_id"].get<G4int>();
        rec.total_photons = 0;

        for (auto& ex : e["excitation_points"]) {
            ExcitationPoint ep;
            // Coordinates in JSON are in cm, relative to hole centre
            ep.x = ex["x"].get<G4double>() * cm;
            ep.y = ex["y"].get<G4double>() * cm;
            ep.z = ex["z"].get<G4double>() * cm;   // z from SiPM (0 = bottom)
            ep.n_photons = ex["n_photons"].get<G4int>();
            rec.total_photons += ep.n_photons;
            rec.excitations.push_back(ep);
        }
        if (rec.total_photons > 0)
            gElectrons.push_back(rec);
    }
    G4cout << "Loaded " << gElectrons.size()
           << " electron records with excitation points." << G4endl;
}

// ─────────────────────────────────────────────────────────────────
//  Detector Construction
// ─────────────────────────────────────────────────────────────────
class ELDetectorConstruction : public G4VUserDetectorConstruction {
public:
    G4LogicalVolume* fHoleLogical  = nullptr;
    G4LogicalVolume* fSiPMLogical  = nullptr;
    G4LogicalVolume* fBlockLogical = nullptr;
    G4LogicalVolume* fWorldLogical = nullptr;

    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist = G4NistManager::Instance();

        // ── Materials ──────────────────────────────────────────
        const G4double photE = (h_Planck * c_light) / (XE_WL_NM * 1e-9 * m);
        std::vector<G4double> eVec  = {photE * 0.99, photE, photE * 1.01};

        // Xenon gas
        G4Material* xenon = nist->FindOrBuildMaterial("G4_Xe");
        G4MaterialPropertiesTable* xeMPT = new G4MaterialPropertiesTable();
        std::vector<G4double> rIdxXe = {1.00069, 1.00069, 1.00069};
        std::vector<G4double> absXe  = {1e5*cm, 1e5*cm, 1e5*cm};
        xeMPT->AddProperty("RINDEX", eVec, rIdxXe);
        xeMPT->AddProperty("ABSLENGTH", eVec, absXe);
        xenon->SetMaterialPropertiesTable(xeMPT);

        // PTFE (hole walls) — NO RINDEX so OpBoundary handles the interface
        G4Material* ptfe = nist->FindOrBuildMaterial("G4_TEFLON");

        // Silicon (SiPM)
        G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
        G4MaterialPropertiesTable* siMPT = new G4MaterialPropertiesTable();
        std::vector<G4double> rIdxSi = {3.5, 3.5, 3.5};
        siMPT->AddProperty("RINDEX", eVec, rIdxSi);
        silicon->SetMaterialPropertiesTable(siMPT);

        // ── World volume (Xe gas) ─────────────────────────────
        G4Box* worldSolid = new G4Box("World", WORLD_HALF, WORLD_HALF, WORLD_HALF);
        fWorldLogical = new G4LogicalVolume(worldSolid, xenon, "World");
        G4VPhysicalVolume* worldPhys = new G4PVPlacement(
            nullptr, G4ThreeVector(), fWorldLogical, "World", nullptr, false, 0);

        // ── PTFE block (solid cylinder, mother of hole + SiPM) ─
        G4double sipm_thick = 0.01 * cm;
        G4double block_half_h = (HOLE_HEIGHT + sipm_thick) / 2;
        G4double block_outer_r = HOLE_RADIUS + WALL_THICK;

        G4Tubs* blockSolid = new G4Tubs("Block",
            0, block_outer_r, block_half_h, 0, 360*deg);
        fBlockLogical = new G4LogicalVolume(blockSolid, ptfe, "Block");
        G4double block_zc = HOLE_HEIGHT / 2 - sipm_thick / 2;
        G4VPhysicalVolume* blockPhys = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, block_zc),
            fBlockLogical, "Block", fWorldLogical, false, 0);

        // ── Hole (Xe cylinder) as daughter of PTFE block ──────
        G4Tubs* holeSolid = new G4Tubs("Hole",
            0, HOLE_RADIUS, HOLE_HEIGHT / 2, 0, 360*deg);
        fHoleLogical = new G4LogicalVolume(holeSolid, xenon, "Hole");
        G4VPhysicalVolume* holePhys = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, sipm_thick / 2),
            fHoleLogical, "Hole", fBlockLogical, false, 0);

        // ── SiPM (Si disk) as daughter of PTFE block ──────────
        G4Tubs* sipmSolid = new G4Tubs("SiPM",
            0, HOLE_RADIUS, sipm_thick / 2, 0, 360*deg);
        fSiPMLogical = new G4LogicalVolume(sipmSolid, silicon, "SiPM");
        G4VPhysicalVolume* sipmPhys = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, -HOLE_HEIGHT / 2),
            fSiPMLogical, "SiPM", fBlockLogical, false, 0);

        // ── PTFE optical surface (diffuse Lambertian reflector) ───
        G4OpticalSurface* ptfeSurf = new G4OpticalSurface("PTFESurface");
        ptfeSurf->SetType(dielectric_metal);
        ptfeSurf->SetModel(glisur);
        ptfeSurf->SetPolish(0.0);  // 0 = fully diffuse (Lambertian)
        G4MaterialPropertiesTable* ptfeSurfMPT = new G4MaterialPropertiesTable();
        std::vector<G4double> reflVec = {PTFE_REFLECT, PTFE_REFLECT, PTFE_REFLECT};
        ptfeSurfMPT->AddProperty("REFLECTIVITY", eVec, reflVec);
        ptfeSurf->SetMaterialPropertiesTable(ptfeSurfMPT);

        // Border surface: Hole → Block (photon exits hole into PTFE wall)
        new G4LogicalBorderSurface("HoleToBlock", holePhys, blockPhys, ptfeSurf);

        // ── SiPM optical surface (perfect absorber) ──────────
        G4OpticalSurface* sipmSurf = new G4OpticalSurface("SiPMSurface");
        sipmSurf->SetType(dielectric_metal);
        sipmSurf->SetModel(glisur);
        sipmSurf->SetPolish(1.0);  // polished
        G4MaterialPropertiesTable* sipmSurfMPT = new G4MaterialPropertiesTable();
        std::vector<G4double> noRefl = {0.0, 0.0, 0.0};
        sipmSurfMPT->AddProperty("REFLECTIVITY", eVec, noRefl);
        sipmSurf->SetMaterialPropertiesTable(sipmSurfMPT);

        // Border surface: Hole → SiPM (photon hits SiPM at bottom)
        new G4LogicalBorderSurface("HoleToSiPM", holePhys, sipmPhys, sipmSurf);

        return worldPhys;
    }
};

// ─────────────────────────────────────────────────────────────────
//  Per-event data
// ─────────────────────────────────────────────────────────────────
struct EventData {
    G4int    electron_id   = -1;
    G4int    photon_id     = -1;
    G4int    n_reflections = 0;
    G4bool   reached_sipm  = false;
    G4bool   absorbed_wall = false;
    G4bool   escaped_top   = false;
    G4double hit_x = 0., hit_y = 0.;
    G4double src_x = 0., src_y = 0., src_z = 0.;
};
EventData gEvData;


// ─────────────────────────────────────────────────────────────────
//  Primary Generator Action
// ─────────────────────────────────────────────────────────────────
class ELPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    G4ParticleGun* fGun;
    size_t fElIdx  = 0;   // electron index
    size_t fExIdx  = 0;   // excitation point index
    size_t fPhIdx  = 0;   // photon index within excitation point
    G4int  fPhotonCount = 0;
public:
    ELPrimaryGeneratorAction() {
        fGun = new G4ParticleGun(1);
        fGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
        const G4double photE = (h_Planck * c_light) / (XE_WL_NM * 1e-9 * m);
        fGun->SetParticleEnergy(photE);
    }
    ~ELPrimaryGeneratorAction() { delete fGun; }

    void GeneratePrimaries(G4Event* event) override {
        if (fElIdx >= gElectrons.size()) return;
        const auto& elec = gElectrons[fElIdx];
        const auto& ep   = elec.excitations[fExIdx];

        // Position at excitation point (already in G4 units)
        fGun->SetParticlePosition(G4ThreeVector(ep.x, ep.y, ep.z));

        // Random isotropic direction
        G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
        G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
        G4ThreeVector dir(sinTheta * std::cos(phi),
                          sinTheta * std::sin(phi),
                          cosTheta);
        fGun->SetParticleMomentumDirection(dir);

        // Polarisation
        G4ThreeVector perp = dir.orthogonal().unit();
        G4double polAngle = G4UniformRand() * 2 * CLHEP::pi;
        G4ThreeVector pol = perp.rotate(polAngle, dir);
        fGun->SetParticlePolarization(pol);

        // Fill event data
        gEvData.electron_id   = elec.id;
        gEvData.photon_id     = fPhotonCount++;
        gEvData.n_reflections = 0;
        gEvData.reached_sipm  = false;
        gEvData.absorbed_wall = false;
        gEvData.escaped_top   = false;
        gEvData.src_x = ep.x / cm;
        gEvData.src_y = ep.y / cm;
        gEvData.src_z = ep.z / cm;
        gEvData.hit_x = 0.;
        gEvData.hit_y = 0.;

        fGun->GeneratePrimaryVertex(event);

        // Advance pointers
        fPhIdx++;
        if ((G4int)fPhIdx >= ep.n_photons) {
            fPhIdx = 0;
            fExIdx++;
            if (fExIdx >= elec.excitations.size()) {
                fExIdx = 0;
                fElIdx++;
            }
        }
    }
    bool Done() const { return fElIdx >= gElectrons.size(); }
};

// ─────────────────────────────────────────────────────────────────
//  Stepping Action
// ─────────────────────────────────────────────────────────────────
class ELSteppingAction : public G4UserSteppingAction {
    G4OpBoundaryProcess* fBoundaryProc = nullptr;

    G4OpBoundaryProcess* FindBoundaryProcess() {
        auto* pm = G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
        auto* pv = pm->GetProcessList();
        for (G4int i = 0; i < (G4int)pv->size(); i++) {
            auto* bp = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
            if (bp) return bp;
        }
        return nullptr;
    }

public:
    void UserSteppingAction(const G4Step* step) override {
        // Lazy-init boundary process pointer
        if (!fBoundaryProc) fBoundaryProc = FindBoundaryProcess();

        const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (!proc) return;
        G4String procName = proc->GetProcessName();

        // Check boundary interactions at geometry boundaries
        // OpBoundary acts AFTER Transportation, so check step status
        G4StepStatus ss = step->GetPostStepPoint()->GetStepStatus();
        if (ss == fGeomBoundary && fBoundaryProc) {
            G4OpBoundaryProcessStatus bStatus = fBoundaryProc->GetStatus();
            if (bStatus == LambertianReflection ||
                bStatus == LobeReflection ||
                bStatus == SpikeReflection ||
                bStatus == BackScattering ||
                bStatus == TotalInternalReflection ||
                bStatus == FresnelReflection) {
                gEvData.n_reflections++;
            } else if (bStatus == Absorption) {
                // Wall absorption (exclude SiPM which also has Absorption status)
                G4VPhysicalVolume* pVol = step->GetPostStepPoint()
                    ->GetTouchableHandle()->GetVolume();
                if (!pVol || pVol->GetName() != "SiPM") {
                    gEvData.absorbed_wall = true;
                }
            }
        }

        // Check if entering SiPM
        G4VPhysicalVolume* postVol = step->GetPostStepPoint()
                                         ->GetTouchableHandle()->GetVolume();
        if (postVol && postVol->GetName() == "SiPM") {
            gEvData.reached_sipm = true;
            G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
            gEvData.hit_x = pos.x() / cm;
            gEvData.hit_y = pos.y() / cm;
        }

        // Bulk absorption
        if (procName == "OpAbsorption") {
            gEvData.absorbed_wall = true;
        }

        // Photon escaped above the hole
        G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();
        if (postPos.z() > HOLE_HEIGHT + 0.1*mm) {
            gEvData.escaped_top = true;
            step->GetTrack()->SetTrackStatus(fStopAndKill);
        }
    }
};

// ─────────────────────────────────────────────────────────────────
//  Event Action — fill ntuple
// ─────────────────────────────────────────────────────────────────
class ELEventAction : public G4UserEventAction {
public:
    void EndOfEventAction(const G4Event*) override {
        auto* am = G4AnalysisManager::Instance();
        am->FillNtupleIColumn(0, gEvData.electron_id);
        am->FillNtupleIColumn(1, gEvData.photon_id);
        am->FillNtupleIColumn(2, gEvData.n_reflections);
        am->FillNtupleIColumn(3, (G4int)gEvData.reached_sipm);
        am->FillNtupleIColumn(4, (G4int)gEvData.absorbed_wall);
        am->FillNtupleIColumn(5, (G4int)gEvData.escaped_top);
        am->FillNtupleDColumn(6, gEvData.src_x);
        am->FillNtupleDColumn(7, gEvData.src_y);
        am->FillNtupleDColumn(8, gEvData.src_z);
        am->FillNtupleDColumn(9, gEvData.hit_x);
        am->FillNtupleDColumn(10, gEvData.hit_y);
        am->AddNtupleRow();
    }
};

// ─────────────────────────────────────────────────────────────────
//  Run Action
// ─────────────────────────────────────────────────────────────────
class ELRunAction : public G4UserRunAction {
    std::string fOutFile;
public:
    ELRunAction(const std::string& outfile) : fOutFile(outfile) {}

    void BeginOfRunAction(const G4Run*) override {
        auto* am = G4AnalysisManager::Instance();
        am->SetVerboseLevel(0);
        am->OpenFile(fOutFile);

        am->CreateNtuple("photons", "EL Photon Tracking in Hole");
        am->CreateNtupleIColumn("electron_id");    // 0
        am->CreateNtupleIColumn("photon_id");      // 1
        am->CreateNtupleIColumn("n_reflections");  // 2
        am->CreateNtupleIColumn("reached_sipm");   // 3
        am->CreateNtupleIColumn("absorbed_wall");  // 4
        am->CreateNtupleIColumn("escaped_top");    // 5
        am->CreateNtupleDColumn("src_x_cm");       // 6
        am->CreateNtupleDColumn("src_y_cm");       // 7
        am->CreateNtupleDColumn("src_z_cm");       // 8
        am->CreateNtupleDColumn("hit_x_cm");       // 9
        am->CreateNtupleDColumn("hit_y_cm");       // 10
        am->FinishNtuple();
    }

    void EndOfRunAction(const G4Run*) override {
        auto* am = G4AnalysisManager::Instance();
        am->Write();
        am->CloseFile();
        G4cout << "\n[ELHolePhotonSim] Results written to: " << fOutFile << G4endl;
    }
};

// ─────────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <el_electrons.json> <output.root>\n";
        return 1;
    }
    const std::string infile  = argv[1];
    const std::string outfile = argv[2];

    try {
        LoadElectronData(infile);
    } catch (const std::exception& ex) {
        std::cerr << "Error loading input: " << ex.what() << "\n";
        return 1;
    }

    // Count total photons
    G4long totalPhotons = 0;
    for (const auto& e : gElectrons)
        totalPhotons += e.total_photons;
    G4cout << "Total photons to simulate: " << totalPhotons << G4endl;

    if (totalPhotons == 0) {
        G4cout << "No photons to simulate. Check input data." << G4endl;
        return 0;
    }

    // ── Geant4 setup ─────────────────────────────────────────────
    G4RunManager* runManager = new G4RunManager();

    // Physics list
    G4VModularPhysicsList* physicsList = new FTFP_BERT();
    G4OpticalPhysics* optPhysics = new G4OpticalPhysics();
    physicsList->RegisterPhysics(optPhysics);

    runManager->SetUserInitialization(new ELDetectorConstruction());
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserAction(new ELPrimaryGeneratorAction());
    runManager->SetUserAction(new ELEventAction());
    runManager->SetUserAction(new ELSteppingAction());
    runManager->SetUserAction(new ELRunAction(outfile));
    runManager->Initialize();

    G4cout << "Starting simulation of " << totalPhotons << " photons..." << G4endl;
    runManager->BeamOn((G4int)totalPhotons);

    delete runManager;
    G4cout << "Simulation complete." << G4endl;
    return 0;
}
