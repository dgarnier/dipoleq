/*
 *  DipoleEq Python bindings — nanobind version
 *
 *  Author: Darren Garnier <darren@openstar.nz>
 *
 *  See PyDipolEq.cpp for the pybind11 equivalent (kept for reference).
 */

#include <cstring>
#include <stdexcept>
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>

#include "nrutil.h"
#include "tokamak.h"
#include "shell.h"
#include "psigrid.h"
#include "CPlasmaModel.h"
#include "plasma.h"
#include "AddCoilJ.h"
#include "AddShellJ.h"
#include "FileInput.h"
#include "LoadBndryGreens.h"
#include "LoadMeasGreens.h"
#include "LeastSquares.h"
#include "FindJ.h"
#include "InitJ.h"
#include "PsiBoundary.h"
#include "PlasmaBoundary.h"
#include "FileOutput.h"
#include "Find_ShellCurrent.h"
#include "FindMeasFit.h"
#include "GetPlasmaParameters.h"
#include "GetFluxMoments.h"
#include "Restart.h"
#include "measurement.h"

#include "NbDipolEq.hpp"
#include "version.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;

// ── nrerror override ─────────────────────────────────────────────────────────
// Throw instead of calling exit() so Python can catch C-level errors.
void nrerror(const char error_text[]) {
    throw std::runtime_error(error_text);
}

extern "C" FILE *LogFile;

// ── Logger ───────────────────────────────────────────────────────────────────
class Logger {
public:
    Logger() {
        m_filename = "/dev/null";
        m_file = fopen(m_filename.c_str(), "w");
        LogFile = m_file;
    }
    void set_filename(const std::string &filename) {
        m_filename = filename;
        if (m_file) fclose(m_file);
        m_file = fopen(m_filename.c_str(), "w");
        LogFile = m_file;
    }
    ~Logger() {
        fclose(m_file);
        m_file = nullptr;
        LogFile = nullptr;
    }
private:
    std::string m_filename;
    FILE *m_file;
};

// ── NbMachine — RAII wrapper around TOKAMAK* ─────────────────────────────────
//
// pybind11 wrapped TOKAMAK directly via a factory + __del__.  In nanobind we
// use a proper C++ RAII class to guarantee free_Tokamak() is always called.
class NbMachine {
public:
    TOKAMAK *tok;

    NbMachine() : tok(new_Tokamak()) {}

    explicit NbMachine(const std::string &fname) {
        tok = FileInput(fname.c_str());
        if (!tok)
            throw std::runtime_error("Failed to load machine from file: " + fname);
    }

    ~NbMachine() {
        if (tok) { free_Tokamak(tok, 0 /*not full*/); tok = nullptr; }
    }

    NbMachine(const NbMachine &) = delete;
    NbMachine &operator=(const NbMachine &) = delete;
    NbMachine(NbMachine &&o) noexcept : tok(o.tok) { o.tok = nullptr; }
};

// ── Helper free-functions (mirrors the ones in PyDipolEq.cpp) ────────────────

static void free_COIL(COIL *c, TOKAMAK *m) {
    free_Coil(c, m->PsiGrid->Nsize);
}

static void free_SHELL(SHELL *c, TOKAMAK *m) {
    free_Shell(c, m->PsiGrid->Nsize, m->NumCoils);
}

static void free_MEAS(MEAS *c, TOKAMAK *m) {
    int nsubshells = 0;
    for (int i = 0; i < m->NumShells; i++)
        if (m->Shells[i]) nsubshells += m->Shells[i]->NumSubShells;
    free_Measure(c, m->PsiGrid->Nsize, m->NumCoils, nsubshells);
}

static void nb_set_NumCoils(NbMachine &self, int n) {
    TOKAMAK &tok = *self.tok;
    if (tok.Coils) {
        for (int i = 0; i < tok.NumCoils; i++)
            if (tok.Coils[i]) free_Coil(tok.Coils[i], tok.PsiGrid->Nsize);
        free(tok.Coils);
    }
    tok.NumCoils = n;
    tok.Coils = (COIL **)malloc(n * sizeof(COIL *));
    for (int i = 0; i < n; i++) tok.Coils[i] = new_Coil(0);
}

static void nb_set_NumShells(NbMachine &self, int n) {
    TOKAMAK &tok = *self.tok;
    if (tok.Shells) {
        for (int i = 0; i < tok.NumShells; i++)
            if (tok.Shells[i]) free_Shell(tok.Shells[i], tok.PsiGrid->Nsize, tok.NumCoils);
        free(tok.Shells);
    }
    tok.NumShells = n;
    tok.Shells = (SHELL **)malloc(n * sizeof(SHELL *));
    for (int i = 0; i < n; i++) tok.Shells[i] = new_Shell(0);
}

static void nb_set_NumSeps(NbMachine &self, int n) {
    TOKAMAK &tok = *self.tok;
    if (tok.Seps) {
        for (int i = 0; i < tok.NumSeps; i++)
            if (tok.Seps[i]) free(tok.Seps[i]);
        free(tok.Seps);
    }
    tok.NumSeps = n;
    tok.Seps = (SEPARATRIX **)malloc(n * sizeof(SEPARATRIX *));
    for (int i = 0; i < n; i++) tok.Seps[i] = new_Separatrix();
}

static void nb_set_NumMeasures(NbMachine &self, int n) {
    TOKAMAK &tok = *self.tok;
    if (tok.Measures) {
        int nsubshells = 0;
        for (int j = 0; j < tok.NumShells; j++)
            if (tok.Shells[j]) nsubshells += tok.Shells[j]->NumSubShells;
        for (int i = 0; i < tok.NumMeasures; i++)
            if (tok.Measures[i])
                free_Measure(tok.Measures[i], tok.PsiGrid->Nsize, tok.NumCoils, nsubshells);
        free(tok.Measures);
    }
    tok.NumMeasures = n;
    tok.Measures = (MEAS **)malloc(n * sizeof(MEAS *));
    for (int i = 0; i < n; i++) tok.Measures[i] = nullptr;
}

static void nb_set_NumLimiters(NbMachine &self, int n) {
    TOKAMAK &tok = *self.tok;
    if (tok.Limiters) {
        for (int i = 0; i < tok.NumLimiters; i++)
            if (tok.Limiters[i]) free_Limiter(tok.Limiters[i]);
        free(tok.Limiters);
    }
    tok.NumLimiters = n;
    tok.Limiters = (LIMITER **)malloc(n * sizeof(LIMITER *));
    for (int i = 0; i < n; i++) tok.Limiters[i] = new_Limiter();
}

static void nb_set_NumSubCoils(COIL &self, int n) {
    if (self.SubCoils) {
        for (int i = 0; i < self.NumSubCoils; i++)
            if (self.SubCoils[i]) free(self.SubCoils[i]);
        free(self.SubCoils);
    }
    self.NumSubCoils = n;
    self.SubCoils = (SUBCOIL **)malloc(n * sizeof(SUBCOIL *));
    for (int i = 0; i < n; i++) self.SubCoils[i] = new_SubCoil();
}

static void nb_set_SHELL_NumSubShells(SHELL &self, int n, TOKAMAK *tok) {
    if (self.SubShells) {
        for (int i = 0; i < self.NumSubShells; i++)
            if (self.SubShells[i])
                free_SubShell(self.SubShells[i], tok->PsiGrid->Nsize, tok->NumCoils);
        free(self.SubShells);
    }
    self.NumSubShells = n;
    self.SubShells = (SUBSHELL **)malloc(n * sizeof(SUBSHELL *));
    for (int i = 0; i < n; i++) self.SubShells[i] = new_SubShell();
}

static void nb_set_TOKAMAK_NumSubShells(NbMachine &self, int i, int n) {
    SHELL *shell = self.tok->Shells[i];
    nb_set_SHELL_NumSubShells(*shell, n, self.tok);
}

// ── GetFluxContour wrapper ────────────────────────────────────────────────────
// Returns a pair of 1-D numpy arrays (r, z) for the contour at normalised psi.
// Uses the same "dirty trick" as PyDipolEq.cpp to stash len before the
// nrutil-allocated vector so the capsule deleter can free it properly.
static std::pair<
    nb::ndarray<nb::numpy, double, nb::ndim<1>>,
    nb::ndarray<nb::numpy, double, nb::ndim<1>>>
nb_get_flux_contour(PSIGRID *pg, double psi) {
    int len;
    double *r, *z;
    GetFluxContour(pg, psi, &r, &z, &len);

    *(int *)(r - 1) = len;
    *(int *)(z - 1) = len;

    auto deleter = [](void *p) noexcept {
        double *v = reinterpret_cast<double *>(p);
        int n = *(int *)(v - 1);
        free_dvector(v, 0, n);
    };

    nb::capsule r_cap(r, deleter);
    nb::capsule z_cap(z, deleter);

    size_t shape[1] = {(size_t)len};
    return {
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(r, 1, shape, r_cap),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(z, 1, shape, z_cap)
    };
}

// ── Module ────────────────────────────────────────────────────────────────────

NB_MODULE(core, m) {
    static Logger logger;
    m.doc() = "Python bindings for DipolEq (nanobind version)";
    m.attr("__version__")      = VERSION_FULL;
    auto version_info = nb::dict();
    version_info["version"]        = VERSION;
    version_info["version_full"]   = VERSION_FULL;
    version_info["git_commit"]     = GIT_COMMIT;
    version_info["git_short_hash"] = GIT_SHORT_HASH;
    version_info["git_distance"]   = GIT_DISTANCE;
    m.attr("__version_info__") = version_info;

    // ── Enums ────────────────────────────────────────────────────────────────
    nb::enum_<ModelType>(m, "ModelType", nb::is_arithmetic(), "Plasma model type")
        .value("Std",              ModelType::Std)
        .value("IsoNoFlow",        ModelType::IsoNoFlow)
        .value("IsoFlow",          ModelType::IsoFlow)
        .value("AnisoNoFlow",      ModelType::AnisoNoFlow)
        .value("AnisoFlow",        ModelType::AnisoFlow)
        .value("DipoleStd",        ModelType::DipoleStd)
        .value("DipoleIntStable",  ModelType::DipoleIntStable)
        .value("DipoleStablePsiN", ModelType::DipoleStablePsiN)
        .export_values();

    nb::enum_<MeasType>(m, "MeasType", nb::is_arithmetic(), "Measurement type")
        .value("unk",       MeasType::unk)
        .value("bp",        MeasType::bp)
        .value("press",     MeasType::press)
        .value("pperp",     MeasType::pperp)
        .value("ppar",      MeasType::ppar)
        .value("flux",      MeasType::flux)
        .value("saddle",    MeasType::saddle)
        .value("circle",    MeasType::circle)
        .value("coilcur",   MeasType::coilcur)
        .value("plasmacur", MeasType::plasmacur)
        .value("bt",        MeasType::bt)
        .value("diam",      MeasType::diam)
        .value("bangle",    MeasType::bangle)
        .value("flowt",     MeasType::flowt)
        .value("flowp",     MeasType::flowp)
        .value("ne",        MeasType::ne)
        .value("Te",        MeasType::Te)
        .value("Ti",        MeasType::Ti)
        .value("Zeff",      MeasType::Zeff)
        .value("rho",       MeasType::rho)
        .value("rot",       MeasType::rot)
        .value("ppsix",     MeasType::ppsix)
        .value("bpangle",   MeasType::bpangle)
        .value("pnorm",     MeasType::pnorm)
        .value("J0",        MeasType::J0)
        .export_values();

    nb::enum_<CircleType>(m, "CircleType", nb::is_arithmetic(), "Circle measurement sub-type")
        .value("btcos", CircleType::btcos)
        .value("brsin", CircleType::brsin)
        .value("brcos", CircleType::brcos)
        .export_values();

    // ── CPlasmaModel ─────────────────────────────────────────────────────────
    nb::class_<CPlasmaModel>(m, "CPlasmaModel", "Plasma model base class")
        .def("update_model", [](CPlasmaModel &self, NbMachine &mach) {
            self.UpdateModel(mach.tok);
        }, "Update the plasma model")
        .def("model_input", [](CPlasmaModel &self,
                                const std::string &key,
                                const std::string &val,
                                const std::string &unit) {
            self.ModelInput(const_cast<char *>(key.c_str()),
                            const_cast<char *>(val.c_str()),
                            const_cast<char *>(unit.c_str()));
        }, "Set a model input parameter")
    ;

    // ── PLASMA ───────────────────────────────────────────────────────────────
    nb::class_<PLASMA>(m, "Plasma", "Plasma data and profile results")
        .def("init",     [](PLASMA &self)           { init_Plasma(&self); },
                         "Initialize Plasma")
        .def("plasmaP",  [](PLASMA &self, double pn){ return PlasmaP(&self, pn); },
                         "Calculate plasma pressure")
        .def("plasmaPp", [](PLASMA &self, double pn){ return PlasmaPp(&self, pn); },
                         "Calculate plasma Pprime")
        .def("plasmaG",  [](PLASMA &self, double pn){ return PlasmaG(&self, pn); },
                         "Calculate plasma G")
        .def("plasmaG2p",[](PLASMA &self, double pn){ return PlasmaG2p(&self, pn); },
                         "Calculate plasma G2prime")

        .def_rw("Nsize",       &PLASMA::Nsize,
                "Size of grid -1, should be factor of 2")
        .def_prop_rw("ModelType",
            [](PLASMA &self)             { return ModelType(self.ModelType); },
            [](PLASMA &self, ModelType t){ self.ModelType = (int)t; },
            "Plasma model type, see ModelType enum")
        .def_prop_ro("Model",
            [](PLASMA &self) -> nb::object {
                if (!self.Model) return nb::none();
                return nb::cast(self.Model, nb::rv_policy::reference);
            }, "The plasma model object")

        .def_rw("R0",           &PLASMA::R0,           "Reference major radius")
        .def_rw("Z0",           &PLASMA::Z0,           "Reference vertical position")
        .def_rw("B0",           &PLASMA::B0,           "Vacuum magnetic field at R0, Z0")
        .def_rw("Ip0",          &PLASMA::Ip0,          "initial plasma current")
        .def_rw("B0R0",         &PLASMA::B0R0,         ":math:`B_0 * R_0`")
        .def_rw("Jedge",        &PLASMA::Jedge,        "Edge current density")
        .def_rw("NumBndMomts",  &PLASMA::NumBndMomts,  "Number of moments of boundary to calculate")
        .def_rw("NumPsiPts",    &PLASMA::NumPsiPts,    "Number of normalized psi points to calculate")
        .def_rw("PsiXmax",      &PLASMA::PsiXmax,      "Outermost normalized Psi from 0.0 to 1.0 (say 0.995)")

        // old plasma terms
        .def_rw("G2pTerms",  &PLASMA::G2pTerms, "Number of G^2' terms in polynomial plasma models")
        .def_rw("HTerms",    &PLASMA::HTerms,   "Number of H terms in polynomial plasma models")
        .def_rw("PpTerms",   &PLASMA::PpTerms,  "Number of p' terms in polynomial plasma models")
        .def_rw("RotTerms",  &PLASMA::RotTerms, "Number of Rot terms in flow polynomial plasma models")
        .def_rw("SisoTerms", &PLASMA::SisoTerms,"Number of Siso terms in isotropic polynomial plasma models")
        .def_rw("SparTerms", &PLASMA::SparTerms,"Number of Spar terms in anisotropic polynomial plasma models")
        .def_rw("SperTerms", &PLASMA::SperTerms,"Number of Sper terms in anisotropic polynomial plasma models")

        .def_prop_ro("G2p",  [](PLASMA &s){ return make_dvector(s.G2pTerms,  s.G2p);  }, "G^2' polynomial terms")
        .def_prop_ro("H",    [](PLASMA &s){ return make_dvector(s.HTerms,    s.H);    }, "H polynomial terms")
        .def_prop_ro("Pp",   [](PLASMA &s){ return make_dvector(s.PpTerms,   s.Pp);   }, "p' polynomial terms")
        .def_prop_ro("Rot",  [](PLASMA &s){ return make_dvector(s.RotTerms,  s.Rot);  }, "Rotational flow terms")
        .def_prop_ro("Siso", [](PLASMA &s){ return make_dvector(s.SisoTerms, s.Siso); }, "Isotropic S polynomial terms")
        .def_prop_ro("Spar", [](PLASMA &s){ return make_dvector(s.SparTerms, s.Spar); }, "Anisotropic S parallel terms")
        .def_prop_ro("Sper", [](PLASMA &s){ return make_dvector(s.SperTerms, s.Sper); }, "Anisotropic S perpendicular terms")

        // results
        .def_ro("Ip",           &PLASMA::Ip,           "Plasma current")
        .def_ro("beta0",        &PLASMA::beta0,        "Vacuum toroidal beta at R0")
        .def_ro("beta",         &PLASMA::beta,         "Average toroidal beta")
        .def_ro("betap",        &PLASMA::betap,        "Poloidal beta")
        .def_ro("li",           &PLASMA::li,           "Normalized internal inductance")
        .def_ro("Ltotal",       &PLASMA::Ltotal,       "Total inductance")
        .def_ro("mu",           &PLASMA::mu,           "Normalized diamagnetism")
        .def_ro("Volume",       &PLASMA::Volume,       "Volume")
        .def_ro("CrossSection", &PLASMA::CrossSection, "Cross section")
        .def_ro("Perimeter",    &PLASMA::Perimeter,    "Perimeter")
        .def_ro("Diamag",       &PLASMA::Diamag,       "Diamagnetism")
        .def_ro("q0",           &PLASMA::q0,           "Central safety factor")
        .def_ro("qCircular",    &PLASMA::qCircular,    "qCircular")
        .def_ro("qStar",        &PLASMA::qStar,        "qStar")
        .def_ro("RMagAxis",     &PLASMA::XMagAxis,     "Magnetic axis R")
        .def_ro("ZMagAxis",     &PLASMA::ZMagAxis,     "Magnetic axis Z")
        .def_ro("PsiMagAxis",   &PLASMA::PsiMagAxis,   "Psi at magnetic axis")
        .def_ro("PsiAxis",      &PLASMA::PsiAxis,      "Psi at axis or inner separatrix")
        .def_ro("PsiLim",       &PLASMA::PsiLim,       "Psi at plasma/vacuum boundary")
        .def_ro("PsiFCFS",      &PLASMA::PsiAxis,      "Psi at axis or inner separatrix")
        .def_ro("PsiLCFS",      &PLASMA::PsiLim,       "Psi at outer separatrix")
        .def_ro("HalfWidth",    &PLASMA::HalfWidth,    "Half width")
        .def_ro("Elongation",   &PLASMA::Elongation,   "Elongation")
        .def_ro("ChiSqr",       &PLASMA::ChiSqr,       "Chi squared")
        .def_ro("totKinEnergy", &PLASMA::TotKinEnergy, "Total kinetic energy")
        .def_ro("totMagEnergy", &PLASMA::TotMagEnergy, "Total magnetic energy")

        // 2-D grid arrays
        .def_prop_ro("B2",       [](PLASMA &s){ return make_dmatrix(s.Nsize, s.B2);       }, "B^2 on grid")
        .def_prop_ro("GradPsiR", [](PLASMA &s){ return make_dmatrix(s.Nsize, s.GradPsiX); }, "GradPsiR(self) -> ndarray\n\n:math:`\\del\\psi/\\delR` on grid\n")
        .def_prop_ro("GradPsiZ", [](PLASMA &s){ return make_dmatrix(s.Nsize, s.GradPsiZ); }, "∂Psi/∂Z on grid")
        .def_prop_ro("GradPsi2", [](PLASMA &s){ return make_dmatrix(s.Nsize, s.GradPsi2); }, ":math:`|d\\psi/dR|^2 + |d\\psi/dZ|^2` on grid")
        .def_prop_ro("Bt",       [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Bt);       }, "B_toroidal on grid")
        .def_prop_ro("G",        [](PLASMA &s){ return make_dmatrix(s.Nsize, s.G);        }, ":math:`G = F/(B_0 R_0)` on grid")
        .def_prop_ro("Rho",      [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Rho);      }, "Rho on grid")
        .def_prop_ro("Piso",     [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Piso);     }, "Isotropic pressure on grid")
        .def_prop_ro("Ppar",     [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Ppar);     }, "Anisotropic p_parallel pressure on grid")
        .def_prop_ro("Pper",     [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Pper);     }, "Anistropic p_perpendicular pressure on grid")
        .def_prop_ro("Alpha",    [](PLASMA &s){ return make_dmatrix(s.Nsize, s.Alpha);    }, "Alpha = :math:`\\mu_0*(p_\\parallel - p_\\perp)/B^2` on grid")

        // 1-D profile arrays
        .def_prop_ro("PsiX_pr",    [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.PsiX_pr);    }, "Normalized Psi ordinate")
        .def_prop_ro("Psi_pr",     [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Psi_pr);     }, "Psi on normalized profile")
        .def_prop_ro("P_pr",       [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.P_pr);       }, "Pressure profile")
        .def_prop_ro("G_pr",       [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.G_pr);       }, "G = (F = B_t R)/(B_0 R_0) profile")
        .def_prop_ro("Pp_pr",      [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Pp_pr);      }, "p' profile")
        .def_prop_ro("G2p_pr",     [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.G2p_pr);     }, "G^2' profile")
        .def_prop_ro("q_pr",       [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.q_pr);       }, "q profile")
        .def_prop_ro("Volp_pr",    [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Volp_pr);    }, "dV/dPsi profile")
        .def_prop_ro("Vol_pr",     [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Vol_pr);     }, "V profile")
        .def_prop_ro("S_pr",       [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.S_pr);       }, "S = global shear of q profile")
        .def_prop_ro("B2_pr",      [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.B2_pr);      }, "<B^2> flux tube averaged profile")
        .def_prop_ro("Well_pr",    [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Well_pr);    }, "Magnetic well profile")
        .def_prop_ro("J_pr",       [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.J_pr);       }, "<J> flux tube averaged current density profile")
        .def_prop_ro("Beta_pr",    [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.Beta_pr);    }, "<β> flux tube averaged profile")
        .def_prop_ro("BetaMax_pr", [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.BetaMax_pr); }, "max(β) on flux tube profile")
        .def_prop_ro("RBetaMax_pr",[](PLASMA &s){ return make_dvector(s.NumPsiPts, s.XBetaMax_pr);}, "R of max(β) profile")
        .def_prop_ro("ZBetaMax_pr",[](PLASMA &s){ return make_dvector(s.NumPsiPts, s.ZBetaMax_pr);}, "Z of max(β) profile")
        .def_prop_ro("BBetaMax_pr",[](PLASMA &s){ return make_dvector(s.NumPsiPts, s.BBetaMax_pr);}, "B of max(β) profile")
        .def_prop_ro("BMax_pr",    [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.BMax_pr);    }, "max(B) on flux tube profile")
        .def_prop_ro("RBMax_pr",   [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.XBMax_pr);   }, "R of max(B) profile")
        .def_prop_ro("ZBMax_pr",   [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.ZBMax_pr);   }, "Z of max(B) profile")
        .def_prop_ro("RRMax_pr",   [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.XXMax_pr);   }, "Max R of flux tube contour profile")
        .def_prop_ro("ZRMax_pr",   [](PLASMA &s){ return make_dvector(s.NumPsiPts, s.ZXMax_pr);   }, "Z at Max R of flux tube contour profile")
    ;

    // ── PSIGRID ───────────────────────────────────────────────────────────────
    nb::class_<PSIGRID>(m, "PsiGrid", "Magnetic flux grid data and PDE solver")
        .def("init",               [](PSIGRID &s){ init_PsiGrid(&s); },
                                   "Initialize PsiGrid")
        .def("go_PDE",             [](PSIGRID &s){ GoPDE(&s); },
                                   "Solve the PDE on the PsiGrid")
        .def("make_psi_symmetric", [](PSIGRID &s){ MakePsiSymmetric(&s); },
                                   "Make Psi symmetric")
        .def("get_new_residual",   [](PSIGRID &s){ GetNewResidual(&s); },
                                   "Get new residual")
        .def("new_solution",       [](PSIGRID &s){ NewSolution(&s); },
                                   "Get new solution")
        .def("new_M_solution",     [](PSIGRID &s){ NewMSolution(&s); },
                                   "Get new M solution")
        .def("get_Psi",
            [](PSIGRID &s, double r, double z){ return GetPsi(&s, r, z); },
            "Get Psi at (R, Z)")
        .def("get_IsPlasma",
            [](PSIGRID &s, double r, double z){ return GetIsPlasma(&s, r, z); },
            "Get IsPlasma indicator at (R, Z)")
        .def("init_J",  [](PSIGRID &s, PLASMA &pl){ InitJ(&s, &pl); },
                        "Initialize J")

        .def_rw("Nsize",          &PSIGRID::Nsize,
                "Size of grid -1, should be factor of 2")
        .def_rw("Symmetric",      &PSIGRID::Symmetric,
                "Is the grid symmetric")
        .def_rw("MaxRes",         &PSIGRID::MaxRes)
        .def_rw("PastMaxRes",     &PSIGRID::PastMaxRes)
        .def_rw("Rmax",           &PSIGRID::Xmax,
                "Maximum R of computational grid")
        .def_rw("Rmin",           &PSIGRID::Xmin,
                "Minimum R of computational grid")
        .def_rw("Zmax",           &PSIGRID::Zmax,
                "Maximum Z of computational grid")
        .def_rw("Zmin",           &PSIGRID::Zmin,
                "Minimum Z of computational grid")
        .def_prop_ro("dr",        [](PSIGRID &s){ return s.dx; },
                     "R grid spacing")
        .def_prop_ro("dz",        [](PSIGRID &s){ return s.dz; },
                     "Z grid spacing")
        .def_rw("BoundError",     &PSIGRID::BoundError)
        .def_rw("BoundThreshold", &PSIGRID::BoundThreshold)
        .def_rw("ResThreshold",   &PSIGRID::ResThreshold)
        .def_rw("UnderRelax1",    &PSIGRID::UnderRelax1)
        .def_rw("UnderRelax2",    &PSIGRID::UnderRelax2)
        .def_ro("PsiAxis",        &PSIGRID::PsiAxis,
                "Psi at FCFS or Magnetic Axis")
        .def_ro("PsiMagAxis",     &PSIGRID::PsiMagAxis,
                "Psi at Magnetic Axis")
        .def_ro("PsiLim",         &PSIGRID::PsiLim,
                "Psi at plasma/vacuum boundary")
        .def_ro("PsiFCFS",        &PSIGRID::PsiAxis,
                "Psi at axis or inner separatrix")
        .def_ro("PsiLCFS",        &PSIGRID::PsiLim,
                "Psi at outer separatrix")
        .def_prop_ro("DelPsi",    [](PSIGRID &s){ return s.DelPsi; },
                     "Change in psi across plasma")
        .def_rw("RMagAxis",       &PSIGRID::XMagAxis,
                "Magnetic axis R")
        .def_rw("ZMagAxis",       &PSIGRID::ZMagAxis,
                "Magnetic axis Z")

        .def_prop_ro("R",        [](PSIGRID &s){ return make_dvector(s.Nsize + 1, s.X);       },
                     "R grid")
        .def_prop_ro("Z",        [](PSIGRID &s){ return make_dvector(s.Nsize + 1, s.Z);       },
                     "Z grid")
        .def_prop_ro("IsPlasma", [](PSIGRID &s){ return make_imatrix(s.Nsize,     s.IsPlasma);},
                     "grid point is considered inside the plasma.")
        .def_prop_ro("Psi",      [](PSIGRID &s){ return make_dmatrix(s.Nsize,     s.Psi);     },
                     "Magnetic flux on grid")
        .def_prop_ro("Current",  [](PSIGRID &s){ return make_dmatrix(s.Nsize,     s.Current); },
                     "Current density on grid * μ_0")
        .def_prop_ro("Residual", [](PSIGRID &s){ return make_dmatrix(s.Nsize,     s.Residual);},
                     "Current density residual on grid")

        .def("get_contour", &nb_get_flux_contour,
            "Get flux contour at normalised psi; returns (r, z) numpy arrays")
    ;

    // ── LIMITER ───────────────────────────────────────────────────────────────
    nb::class_<LIMITER>(m, "Limiter", "Limiter line-segment definition")
        .def_rw("R1",      &LIMITER::X1)
        .def_rw("Z1",      &LIMITER::Z1)
        .def_rw("R2",      &LIMITER::X2)
        .def_rw("Z2",      &LIMITER::Z2)
        .def_rw("Enabled", &LIMITER::Enabled,
                "1=outer limiter, 0=disabled, -1=inner limiter")
        .def_prop_ro("PsiLim", [](LIMITER &l){ return l.PsiMin; },
                     "Psi at limiting point")
        .def_prop_ro("RLim",   [](LIMITER &l){ return l.Xmin;   },
                     "R at limiting point")
        .def_prop_ro("ZLim",   [](LIMITER &l){ return l.Zmin;   },
                     "Z at limiting point")
        .def_prop_rw("Name",
            [](LIMITER &s)                       { return std::string(s.Name); },
            [](LIMITER &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(LIMITER::Name) - 1);
                s.Name[sizeof(LIMITER::Name) - 1] = '\0';
            }, "Name of the limiter")
    ;

    // ── SUBCOIL ───────────────────────────────────────────────────────────────
    nb::class_<SUBCOIL>(m, "SubCoil", "Individual coil element")
        .def_rw("R",        &SUBCOIL::X)
        .def_rw("Z",        &SUBCOIL::Z)
        .def_rw("Fraction", &SUBCOIL::CurrentFraction, "Fraction of current")
        .def_prop_rw("Name",
            [](SUBCOIL &s)                       { return std::string(s.Name); },
            [](SUBCOIL &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(SUBCOIL::Name) - 1);
                s.Name[sizeof(SUBCOIL::Name) - 1] = '\0';
            }, "Name of the subcoil")
    ;

    nb::class_<ObjVecView<SUBCOIL>>(m, "SubCoils", "Sequence view of SubCoil objects")
        .def("__len__",     [](ObjVecView<SUBCOIL> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<SUBCOIL> &s, size_t i) -> SUBCOIL * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__",
            [](ObjVecView<SUBCOIL> &s, size_t i, SUBCOIL *c){ s[i] = c; })
        .def("new_subcoil",
            [](ObjVecView<SUBCOIL> &) -> SUBCOIL * { return new_SubCoil(); },
            nb::rv_policy::reference, "Create a new SubCoil (must be inserted into Coil)")
    ;

    // ── COIL ─────────────────────────────────────────────────────────────────
    nb::class_<COIL>(m, "Coil", "Coil data and sub-coil management")
        .def_rw("Enabled",     &COIL::Enabled)
        .def_rw("CoilCurrent", &COIL::CoilCurrent)
        .def_rw("R",           &COIL::X,  "Coil centroid R")
        .def_rw("dR",          &COIL::dX, "Coil radial width, > 0 for automatic subcoil generation")
        .def_rw("Z",           &COIL::Z,  "Coil centroid Z")
        .def_rw("dZ",          &COIL::dZ, "Coil vertical width")
        .def_prop_rw("Name",
            [](COIL &s)                       { return std::string(s.Name); },
            [](COIL &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(COIL::Name) - 1);
                s.Name[sizeof(COIL::Name) - 1] = '\0';
            }, "Name of the coil")
        .def_prop_ro("SubCoils",
            [](COIL &s) -> nb::object {
                if (!s.SubCoils) return nb::none();
                return nb::cast(new ObjVecView<SUBCOIL>(s.NumSubCoils, s.SubCoils),
                                nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Return vector of SubCoils")
        .def_prop_rw("NumSubCoils",
            [](COIL &s){ return s.NumSubCoils; },
            &nb_set_NumSubCoils,
            "Number of subcoils, setting will erase old subcoils")
        .def("compute_SubCoils",
            [](COIL &s, PSIGRID &pg){ compute_SubCoils(&s, &pg); },
            "Generate subcoils from coil parameters and grid")
    ;

    // ── SUBSHELL ──────────────────────────────────────────────────────────────
    nb::class_<SUBSHELL>(m, "SubShell", "Conducting shell segment")
        .def_rw("R",       &SUBSHELL::X)
        .def_rw("Z",       &SUBSHELL::Z)
        .def_rw("Radius",  &SUBSHELL::Radius)
        .def_rw("Current", &SUBSHELL::Current)
        .def_prop_rw("Name",
            [](SUBSHELL &s)                       { return std::string(s.Name); },
            [](SUBSHELL &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(SUBSHELL::Name) - 1);
                s.Name[sizeof(SUBSHELL::Name) - 1] = '\0';
            }, "Name of the subshell")
    ;

    nb::class_<ObjVecView<SUBSHELL>>(m, "SubShells", "Sequence view of SubShell objects")
        .def("__len__",    [](ObjVecView<SUBSHELL> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<SUBSHELL> &s, size_t i) -> SUBSHELL * { return s[i]; },
            nb::rv_policy::reference)
    ;

    // ── SHELL ─────────────────────────────────────────────────────────────────
    nb::class_<SHELL>(m, "Shell", "Conducting shell data")
        .def_rw("Enabled",       &SHELL::Enabled)
        .def_ro("NumSubShells",  &SHELL::NumSubShells)
        .def_prop_rw("Name",
            [](SHELL &s)                       { return std::string(s.Name); },
            [](SHELL &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(SHELL::Name) - 1);
                s.Name[sizeof(SHELL::Name) - 1] = '\0';
            })
        .def_prop_ro("SubShells",
            [](SHELL &s) -> nb::object {
                if (!s.SubShells) return nb::none();
                return nb::cast(new ObjVecView<SUBSHELL>(s.NumSubShells, s.SubShells),
                                nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Return vector of SubShells")
        .def("set_NumSubShells",
            [](SHELL &s, int n, NbMachine &mach){
                nb_set_SHELL_NumSubShells(s, n, mach.tok);
            }, "Reallocate sub-shells (requires parent machine for Nsize/NumCoils)")
    ;

    // ── SEPARATRIX ────────────────────────────────────────────────────────────
    nb::class_<SEPARATRIX>(m, "Separatrix", "Separatrix search definition")
        .def_rw("Enabled",      &SEPARATRIX::Enabled)
        .def_rw("IsSeparatrix", &SEPARATRIX::IsSeparatrix)
        .def_rw("R1",           &SEPARATRIX::X1, "Box to search for separatrix")
        .def_rw("Z1",           &SEPARATRIX::Z1, "Box to search for separatrix")
        .def_rw("R2",           &SEPARATRIX::X2, "Box to search for separatrix")
        .def_rw("Z2",           &SEPARATRIX::Z2, "Box to search for separatrix")
        .def_rw("RC",           &SEPARATRIX::XC, "Center of plasma from Sep")
        .def_rw("ZC",           &SEPARATRIX::ZC, "Center of plasma from Sep")
        .def_ro("Psi",          &SEPARATRIX::PsiSep, "Value of Psi at separatrix")
        .def_ro("Rs",           &SEPARATRIX::Xs, "R of separatrix")
        .def_ro("Zs",           &SEPARATRIX::Zs, "Z of separatrix")
        .def_prop_rw("Name",
            [](SEPARATRIX &s)                       { return std::string(s.Name); },
            [](SEPARATRIX &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(SEPARATRIX::Name) - 1);
                s.Name[sizeof(SEPARATRIX::Name) - 1] = '\0';
            })
    ;

    // ── MEAS ─────────────────────────────────────────────────────────────────
    nb::class_<MEAS>(m, "Measure", "Diagnostic measurement data")
        .def_rw("R",      &MEAS::X)
        .def_rw("Z",      &MEAS::Z)
        .def_rw("Value",  &MEAS::Value)
        .def_rw("StdDev", &MEAS::StdDev)
        .def_rw("Fit",    &MEAS::Fit)
        .def_rw("Now",    &MEAS::Now)
        .def_prop_rw("Type",
            [](MEAS &s)          { return MeasType(s.mType); },
            [](MEAS &s, MeasType t){ s.mType = (int)t; },
            "Type of the measurement")
        .def_prop_rw("Name",
            [](MEAS &s)                       { return std::string(s.Name); },
            [](MEAS &s, const std::string &n) {
                strncpy(s.Name, n.c_str(), sizeof(MEAS::Name) - 1);
                s.Name[sizeof(MEAS::Name) - 1] = '\0';
            }, "Name of the measurement")
        // Union-field accessors — same as pybind11 version
        .def_prop_rw("Angle",
            [](MEAS &s)          { return s.parm.bp.Angle; },
            [](MEAS &s, double v){ s.parm.bp.Angle = v; },
            "Angle")
        .def_prop_rw("Radius",
            [](MEAS &s)          { return s.parm.circle.Radius; },
            [](MEAS &s, double v){ s.parm.circle.Radius = v; },
            "Radius")
        .def_prop_rw("Number",
            [](MEAS &s)       { return s.parm.circle.Number; },
            [](MEAS &s, int v){ s.parm.circle.Number = v; },
            "Number")
        .def_prop_rw("CircleType",
            [](MEAS &s)              { return CircleType(s.parm.circle.CircleType); },
            [](MEAS &s, CircleType t){ s.parm.circle.CircleType = (int)t; },
            "CircleType")
        .def_prop_rw("CoilNum",
            [](MEAS &s)       { return s.parm.coilcur.CoilNum; },
            [](MEAS &s, int v){ s.parm.coilcur.CoilNum = v; },
            "CoilNum")
        .def_prop_rw("R1",
            [](MEAS &s)          { return s.parm.saddle.X1; },
            [](MEAS &s, double v){ s.parm.saddle.X1 = v; },
            "Saddle Loop R1")
        .def_prop_rw("Z1",
            [](MEAS &s)          { return s.parm.saddle.Z1; },
            [](MEAS &s, double v){ s.parm.saddle.Z1 = v; },
            "Saddle Loop Z1")
        .def_prop_rw("R2",
            [](MEAS &s)          { return s.parm.saddle.X2; },
            [](MEAS &s, double v){ s.parm.saddle.X2 = v; },
            "Saddle Loop R2")
        .def_prop_rw("Z2",
            [](MEAS &s)          { return s.parm.saddle.Z2; },
            [](MEAS &s, double v){ s.parm.saddle.Z2 = v; },
            "Saddle Loop Z2")
    ;

    // ── ObjVecView specialisations ────────────────────────────────────────────

    nb::class_<ObjVecView<COIL>>(m, "Coils", "Sequence view of Coil objects")
        .def("__len__",    [](ObjVecView<COIL> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<COIL> &s, size_t i) -> COIL * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__", &ObjVecView<COIL>::insert)
        .def("new_Coil",
            [](ObjVecView<COIL> &, int n) -> COIL * { return new_Coil(n); },
            nb::rv_policy::reference, "Create a new Coil (must be inserted via __setitem__)")
    ;

    nb::class_<ObjVecView<LIMITER>>(m, "Limiters", "Sequence view of Limiter objects")
        .def("__len__",    [](ObjVecView<LIMITER> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<LIMITER> &s, size_t i) -> LIMITER * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__",
            [](ObjVecView<LIMITER> &s, size_t i, LIMITER *l){ s[i] = l; })
        .def("new_limiter",
            [](ObjVecView<LIMITER> &) -> LIMITER * { return new_Limiter(); },
            nb::rv_policy::reference, "Add a new limiter")
    ;

    nb::class_<ObjVecView<SEPARATRIX>>(m, "Separatrices", "Sequence view of Separatrix objects")
        .def("__len__",    [](ObjVecView<SEPARATRIX> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<SEPARATRIX> &s, size_t i) -> SEPARATRIX * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__",
            [](ObjVecView<SEPARATRIX> &s, size_t i, SEPARATRIX *p){ s[i] = p; })
        .def("new_separatrix",
            [](ObjVecView<SEPARATRIX> &) -> SEPARATRIX * { return new_Separatrix(); },
            nb::rv_policy::reference, "Add a new separatrix")
    ;

    nb::class_<ObjVecView<MEAS>>(m, "Measures", "Sequence view of Measure objects")
        .def("__len__",    [](ObjVecView<MEAS> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<MEAS> &s, size_t i) -> MEAS * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__",
            [](ObjVecView<MEAS> &s, size_t i, MEAS *p){ s[i] = p; })
        .def("new_meas",
            [](ObjVecView<MEAS> &, MeasType t) -> MEAS * {
                return new_Measure((int)t);
            }, nb::rv_policy::reference,
            "Create a new Measure of the given MeasType")
    ;

    nb::class_<ObjVecView<SHELL>>(m, "Shells", "Sequence view of Shell objects")
        .def("__len__",    [](ObjVecView<SHELL> &s){ return s.size(); })
        .def("__getitem__",
            [](ObjVecView<SHELL> &s, size_t i) -> SHELL * { return s[i]; },
            nb::rv_policy::reference)
        .def("__setitem__", &ObjVecView<SHELL>::insert)
        .def("new_shell",
            [](ObjVecView<SHELL> &) -> SHELL * { return new_Shell(0); },
            nb::rv_policy::reference, "Add a new shell")
    ;

    // ── Machine (NbMachine) ───────────────────────────────────────────────────
    nb::class_<NbMachine>(m, "Machine",
        "Top-level DipolEq machine object with proper RAII lifetime management")
        // constructors
        .def(nb::init<>(), "Create an empty Machine")
        .def(nb::init<std::string>(), "Load a Machine from an input file")

        // initialisation
        .def("init",
            [](NbMachine &s){ init_Tokamak(s.tok); },
            "Initialize Machine")
        .def("set_coil_NumSubCoils",
            [](NbMachine &s, int i, int n){
                if (s.tok->Coils[i]) free_Coil(s.tok->Coils[i], s.tok->PsiGrid->Nsize);
                s.tok->Coils[i] = new_Coil(n);
            }, "Replace coil[i] with a new Coil having n sub-coils")
        .def("set_NumSubcoils",   // alias kept for compatibility
            [](NbMachine &s, int i, int n){
                if (s.tok->Coils[i]) free_Coil(s.tok->Coils[i], s.tok->PsiGrid->Nsize);
                s.tok->Coils[i] = new_Coil(n);
            })
        .def("set_NumSubShells",  &nb_set_TOKAMAK_NumSubShells,
            "Reallocate sub-shells for Shells[i]")

        // timing
        .def("set_start_time", [](NbMachine &s){ SetStartTime(s.tok); },
             "Set the start time")
        .def("set_stop_time",  [](NbMachine &s){ SetStopTime(s.tok); },
             "Set the end time")

        // restart
        .def("read_restart",   [](NbMachine &s){ ReadRestart(s.tok->RSname, s.tok); },
             "Read the restart file")
        .def("write_restart",  [](NbMachine &s){ WriteRestart(s.tok->RSname, s.tok); },
             "Write the restart file")

        // physics operations
        .def("find_shell_current",   [](NbMachine &s){ Find_ShellCurrent(s.tok); },
             "Find shell current")
        .def("load_bndry_greens",    [](NbMachine &s){ LoadBndryGreens(s.tok); },
             "Load boundary greens")
        .def("free_bndry_greens",    [](NbMachine &s){ free_BndryGreens(s.tok); },
             "Free boundary greens")
        .def("psi_boundary",         [](NbMachine &s){ PsiBoundary(s.tok); },
             "Calculate psi boundary")
        .def("add_coil_J",           [](NbMachine &s){ AddCoilJ(s.tok); },
             "Add Coil currents")
        .def("add_shell_J",          [](NbMachine &s){ AddShellJ(s.tok); },
             "Add Shell currents")
        .def("find_plasma_boundary", [](NbMachine &s){ PlasmaBoundary(s.tok); },
             "Find plasma boundary")
        .def("load_meas_greens",     [](NbMachine &s){ LoadMeasGreens(s.tok); },
             "Load measurement greens")
        .def("free_meas_greens",     [](NbMachine &s){ free_MeasGreens(s.tok); },
             "Free measurement greens")
        .def("least_squares",        [](NbMachine &s, int isFirstTime){ LeastSquares(s.tok, isFirstTime); },
             "Least squares")
        .def("get_plasma_parameters",[](NbMachine &s){ GetPlasmaParameters(s.tok); },
             "Get plasma parameters")
        .def("zero_J",               [](NbMachine &s){ ZeroJ(s.tok); },
             "Zero J")
        .def("find_J",               [](NbMachine &s){ FindJ(s.tok); },
             "Find J")
        .def("write_GS2_geo",        [](NbMachine &s){ GS2Output(s.tok); },
             "Write GS2 geometry file")

        // iteration counters / flags
        .def_prop_rw("MaxIterFixed",   [](NbMachine &s){ return s.tok->MaxIterFixed; },
                                       [](NbMachine &s, int v){ s.tok->MaxIterFixed = v; })
        .def_prop_rw("MaxIterFree",    [](NbMachine &s){ return s.tok->MaxIterFree; },
                                       [](NbMachine &s, int v){ s.tok->MaxIterFree = v; })
        .def_prop_rw("IterFixed",      [](NbMachine &s){ return s.tok->IterFixed; },
                                       [](NbMachine &s, int v){ s.tok->IterFixed = v; })
        .def_prop_rw("IterFree",       [](NbMachine &s){ return s.tok->IterFree; },
                                       [](NbMachine &s, int v){ s.tok->IterFree = v; })
        .def_prop_rw("LHGreenStatus",  [](NbMachine &s){ return s.tok->LHGreenStatus; },
                                       [](NbMachine &s, int v){ s.tok->LHGreenStatus = v; })
        .def_prop_rw("MGreenStatus",   [](NbMachine &s){ return s.tok->MGreenStatus; },
                                       [](NbMachine &s, int v){ s.tok->MGreenStatus = v; })
        .def_prop_rw("SGreenStatus",   [](NbMachine &s){ return s.tok->SGreenStatus; },
                                       [](NbMachine &s, int v){ s.tok->SGreenStatus = v; })
        .def_prop_rw("SInductStatus",  [](NbMachine &s){ return s.tok->SInductStatus; },
                                       [](NbMachine &s, int v){ s.tok->SInductStatus = v; })
        .def_prop_rw("RestartStatus",  [](NbMachine &s){ return s.tok->RestartStatus; },
                                       [](NbMachine &s, int v){ s.tok->RestartStatus = v; })
        .def_prop_rw("RestartUnkns",   [](NbMachine &s){ return s.tok->RestartUnkns; },
                                       [](NbMachine &s, int v){ s.tok->RestartUnkns = v; })
        .def_prop_rw("VacuumOnly",     [](NbMachine &s){ return s.tok->VacuumOnly; },
                                       [](NbMachine &s, int v){ s.tok->VacuumOnly = v; })
        .def_prop_rw("NumEqualEq",     [](NbMachine &s){ return s.tok->NumEqualEq; },
                                       [](NbMachine &s, int v){ s.tok->NumEqualEq = v; })
        .def_prop_rw("NumMCarloEq",    [](NbMachine &s){ return s.tok->NumMCarloEq; },
                                       [](NbMachine &s, int v){ s.tok->NumMCarloEq = v; })
        .def_prop_rw("NumMCarloData",  [](NbMachine &s){ return s.tok->NumMCarloData; },
                                       [](NbMachine &s, int v){ s.tok->NumMCarloData = v; })
        .def_prop_rw("MaxIterMCarlo",  [](NbMachine &s){ return s.tok->MaxIterMCarlo; },
                                       [](NbMachine &s, int v){ s.tok->MaxIterMCarlo = v; })
        .def_prop_rw("Confidence",     [](NbMachine &s){ return s.tok->Confidence; },
                                       [](NbMachine &s, double v){ s.tok->Confidence = v; })

        // read-only timestamp strings
        .def_prop_ro("Start", [](NbMachine &s){ return std::string(s.tok->Start); })
        .def_prop_ro("Stop",  [](NbMachine &s){ return std::string(s.tok->Stop);  })

        // sub-objects (non-owning references; lifetime tied to Machine)
        .def_prop_ro("PsiGrid",
            [](NbMachine &s) -> PSIGRID * { return s.tok->PsiGrid; },
            nb::rv_policy::reference_internal, "The PsiGrid object")
        .def_prop_ro("Plasma",
            [](NbMachine &s) -> PLASMA * { return s.tok->Plasma; },
            nb::rv_policy::reference_internal, "The Plasma object")

        // string fields (fixed char arrays in TOKAMAK)
        .def_prop_rw("Name",
            [](NbMachine &s){ return std::string(s.tok->Name); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->Name, v.c_str(), sizeof(TOKAMAK::Name) - 1); },
            "Set the name")
        .def_prop_rw("Info",
            [](NbMachine &s){ return std::string(s.tok->Info); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->Info, v.c_str(), sizeof(TOKAMAK::Info) - 1); },
            "Info about the machine")
        .def_prop_rw("Iname",
            [](NbMachine &s){ return std::string(s.tok->Iname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->Iname, v.c_str(), sizeof(TOKAMAK::Iname) - 1); },
            "Input name")
        .def_prop_rw("Oname",
            [](NbMachine &s){ return std::string(s.tok->Oname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->Oname, v.c_str(), sizeof(TOKAMAK::Oname) - 1); },
            "Output name")
        .def_prop_rw("LHname",
            [](NbMachine &s){ return std::string(s.tok->LHname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->LHname, v.c_str(), sizeof(TOKAMAK::LHname) - 1); },
            "LH name")
        .def_prop_rw("MGname",
            [](NbMachine &s){ return std::string(s.tok->MGname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->MGname, v.c_str(), sizeof(TOKAMAK::MGname) - 1); },
            "MG name")
        .def_prop_rw("SGname",
            [](NbMachine &s){ return std::string(s.tok->SGname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->SGname, v.c_str(), sizeof(TOKAMAK::SGname) - 1); },
            "SG name")
        .def_prop_rw("SMname",
            [](NbMachine &s){ return std::string(s.tok->SMname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->SMname, v.c_str(), sizeof(TOKAMAK::SMname) - 1); },
            "SM name")
        .def_prop_rw("RSname",
            [](NbMachine &s){ return std::string(s.tok->RSname); },
            [](NbMachine &s, const std::string &v){
                strncpy(s.tok->RSname, v.c_str(), sizeof(TOKAMAK::RSname) - 1); },
            "RS name")

        // collection accessors
        .def_prop_ro("Coils",
            [](NbMachine &s) -> nb::object {
                if (!s.tok->Coils) return nb::none();
                return nb::cast(
                    new ObjVecView<COIL>(s.tok->NumCoils, s.tok->Coils, s.tok, free_COIL),
                    nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Array of COILS")
        .def_prop_rw("NumCoils",
            [](NbMachine &s){ return s.tok->NumCoils; },
            &nb_set_NumCoils,
            "Number of coils, setting will erase old coils")

        .def_prop_ro("Limiters",
            [](NbMachine &s) -> nb::object {
                if (!s.tok->Limiters) return nb::none();
                return nb::cast(
                    new ObjVecView<LIMITER>(s.tok->NumLimiters, s.tok->Limiters),
                    nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Get the limiters")
        .def_prop_rw("NumLimiters",
            [](NbMachine &s){ return s.tok->NumLimiters; },
            &nb_set_NumLimiters,
            "Number of limiters, setting will erase old limiters")

        .def_prop_ro("Shells",
            [](NbMachine &s) -> nb::object {
                if (!s.tok->Shells) return nb::none();
                return nb::cast(
                    new ObjVecView<SHELL>(s.tok->NumShells, s.tok->Shells, s.tok, free_SHELL),
                    nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Get the shells")
        .def_prop_rw("NumShells",
            [](NbMachine &s){ return s.tok->NumShells; },
            &nb_set_NumShells,
            "Number of shells, setting will erase old shells")

        .def_prop_ro("Measures",
            [](NbMachine &s) -> nb::object {
                if (!s.tok->Measures) return nb::none();
                return nb::cast(
                    new ObjVecView<MEAS>(s.tok->NumMeasures, s.tok->Measures, s.tok, free_MEAS),
                    nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Get the measurements")
        .def_prop_rw("NumMeasures",
            [](NbMachine &s){ return s.tok->NumMeasures; },
            &nb_set_NumMeasures,
            "Number of measurements, setting will erase old measurements")

        .def_prop_ro("Seps",
            [](NbMachine &s) -> nb::object {
                if (!s.tok->Seps) return nb::none();
                return nb::cast(
                    new ObjVecView<SEPARATRIX>(s.tok->NumSeps, s.tok->Seps),
                    nb::rv_policy::take_ownership);
            }, nb::keep_alive<0, 1>(), "Get the separatrixes")
        .def_prop_rw("NumSeps",
            [](NbMachine &s){ return s.tok->NumSeps; },
            &nb_set_NumSeps,
            "Number of separatrixes, setting will erase old separatrixes")
    ;
}
