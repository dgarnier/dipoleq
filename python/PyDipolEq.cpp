/*
 *  DipoleEq python bindings
 *
 *  Author: Darren Garnier <darren@openstar.nz>
 *
 *
 */


#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

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
#include "Find_ShellCurrent.h"
#include "FindMeasFit.h"
#include "GetPlasmaParameters.h"
#include "GetFluxMoments.h"
#include "Restart.h"
#include "measurement.h"

#include "FileInput.h"
#include "PyDipolEq.hpp"
#include "version.hpp"


namespace py = pybind11;
using namespace pybind11::literals;

// make a minimalistic "C" type binding for dipoleq
// not very pythonic, and not very pybind11.
// TOKAMAK, PSIGIRD, PLASMA, SHELL, COIL, MEASUREMENT
// are all C structs.

// instead of exiting, throw an exception for python to catch
// not very pythonic, but it works.
void nrerror(const char error_text[])
{
	// nrinfo(error_text);  // maybe not even do this
	throw std::runtime_error(error_text);
}

extern "C" FILE *LogFile;

class Logger {
public:
    Logger() {
        m_filename = "/dev/null";
        m_file = fopen(m_filename.c_str(), "w");
        LogFile = m_file;
    };
    void set_filename(std::string filename) {
        m_filename = filename;
        if (m_file) {
            fclose(m_file);
        }
        m_file = fopen(m_filename.c_str(), "w");
        LogFile = m_file;
    };
    ~Logger() {
        fclose(m_file);
        m_file = NULL;
        LogFile = NULL;
    };

private:
    std::string m_filename;
    FILE *m_file;
};


// because nrutils suck and you can't just free a thing
// really this all needs to be rewritten in C++
// with std::vector and eigen

void free_COIL(COIL * c, TOKAMAK *m) {
    int nmax = m->PsiGrid->Nsize;
    free_Coil(c, nmax);
}
// don't need special free_SUBCOIL
// don't need special free_LIMITER
// don't need special free_SEPARATRIX

void free_SHELL(SHELL * c, TOKAMAK *m) {
    int nmax = m->PsiGrid->Nsize;
    int ncoils = m->NumCoils;
    free_Shell(c, nmax, ncoils);
}

// can't actually use this function.
void free_SUBSHELL(SUBSHELL * c, TOKAMAK *m) {
    int nmax = m->PsiGrid->Nsize;
    int ncoils = m->NumCoils;
    free_SubShell(c, nmax, ncoils);
}

void free_MEAS(MEAS * c, TOKAMAK *m) {
    int nmax = m->PsiGrid->Nsize;
    int ncoils = m->NumCoils;
    int nsubshells = 0;
    for (int i=0; i<m->NumShells; i++)
        if (m->Shells[i] != NULL)
            nsubshells += m->Shells[i]->NumSubShells;
    free_Measure(c, nmax, ncoils, nsubshells);
}

void set_NumCoils(TOKAMAK& self, int n) {
    if (self.Coils != NULL) {
        for (int i=0; i<self.NumCoils; i++)
            if (self.Coils[i] != NULL)
                free_Coil(self.Coils[i], self.PsiGrid->Nsize);
        free(self.Coils);
    }
    self.NumCoils = n;
    self.Coils = (COIL **)malloc(n*sizeof(COIL *));
    for (int i=0; i<n; i++)
        self.Coils[i] = new_Coil(0);
}

void set_NumShells(TOKAMAK& self, int n) {
    if (self.Shells != NULL) {
        for (int i=0; i<self.NumShells; i++)
            if (self.Shells[i] != NULL)
                free_Shell(self.Shells[i], self.PsiGrid->Nsize, self.NumCoils);
        free(self.Shells);
    }
    self.NumShells = n;
    self.Shells = (SHELL **)malloc(n*sizeof(SHELL *));
    for (int i=0; i<n; i++)
        self.Shells[i] = new_Shell(0);
}

void set_NumSeps(TOKAMAK& self, int n) {
    if (self.Seps != NULL) {
        for (int i=0; i<self.NumSeps; i++)
            if (self.Seps[i] != NULL)
                free(self.Seps[i]);
        free(self.Seps);
    }
    self.NumSeps = n;
    self.Seps = (SEPARATRIX **)malloc(n*sizeof(SEPARATRIX *));
    for (int i=0; i<n; i++)
        self.Seps[i] = new_Separatrix();
}

void set_NumMeasures(TOKAMAK& self, int n) {
    if (self.Measures != NULL) {
        for (int i=0; i<self.NumMeasures; i++)
            if (self.Measures[i] != NULL) {
                int nsubshells = 0;
                for (int j=0; j<self.NumShells; j++)
                    if (self.Shells[j] != NULL)
                        nsubshells += self.Shells[j]->NumSubShells;
                free_Measure(self.Measures[i],self.PsiGrid->Nsize, self.NumCoils, nsubshells);
            }
        free(self.Measures);
    }
    self.NumMeasures = n;
    self.Measures = (MEAS **)malloc(n*sizeof(MEAS *));
    for (int i=0; i<n; i++)
        self.Measures[i] = NULL;
}

void set_NumLimiters(TOKAMAK& self, int n) {
    if (self.Limiters != NULL) {
        for (int i=0; i<self.NumLimiters; i++)
            if (self.Limiters[i] != NULL)
                free_Limiter(self.Limiters[i]);
        free(self.Limiters);
    }
    self.NumLimiters = n;
    self.Limiters = (LIMITER **)malloc(n*sizeof(LIMITER *));
    for (int i=0; i<n; i++)
        self.Limiters[i] = new_Limiter();
}

void set_TOKAMAK_NumSubShells(TOKAMAK& self, int i, int n) {
    SHELL *shell = self.Shells[i];
    if (shell->SubShells != NULL) {
        for (int i=0; i<shell->NumSubShells; i++)
            if (shell->SubShells[i] != NULL)
                free_SubShell(shell->SubShells[i], self.PsiGrid->Nsize, self.NumCoils);
        free(shell->SubShells);
    }
    shell->NumSubShells = n;
    shell->SubShells = (SUBSHELL **)malloc(n*sizeof(SUBSHELL *));
    for (int i=0; i<n; i++)
        shell->SubShells[i] = new_SubShell();
}

void set_SHELL_NumSubShells(SHELL& self, int n, TOKAMAK *tok) {
    if (self.SubShells != NULL) {
        for (int i=0; i<self.NumSubShells; i++)
            if (self.SubShells[i] != NULL)
                free_SubShell(self.SubShells[i], tok->PsiGrid->Nsize, tok->NumCoils);
        free(self.SubShells);
    }
    self.NumSubShells = n;
    self.SubShells = (SUBSHELL **)malloc(n*sizeof(SUBSHELL *));
    for (int i=0; i<n; i++)
        self.SubShells[i] = new_SubShell();
}

void set_NumSubCoils(COIL& self, int n) {
    if (self.SubCoils != NULL) {
        for (int i=0; i<self.NumSubCoils; i++)
            if (self.SubCoils[i] != NULL)
                free(self.SubCoils[i]);
        free(self.SubCoils);
    }
    self.NumSubCoils = n;
    self.SubCoils = (SUBCOIL **)malloc(n*sizeof(SUBCOIL *));
    for (int i=0; i<n; i++)
        self.SubCoils[i] = new_SubCoil();
}

std::pair<py::array_t<double>, py::array_t<double>> 
get_flux_contour(PSIGRID *pg, double psi) {
    int len;
    double *r, *z;
    // internally calls dvector(0, len) for r and z
    // these suck to get rid of.. I'll cheat and put
    // the length in the NR_END area, which is the
    // root of the problem anyway.
    // FIXME
    // should really fix nrutils to use the area
    // to clean up after themselves and not require 
    // sizes on their destructors.
    GetFluxContour(pg, psi, &r, &z, &len);
    *(int *)(r-1)=len; // oh so dirty
    *(int *)(z-1)=len;

    void (*dvlen)(void *f) = [](void *f) {
        double *v = reinterpret_cast<double *>(f);
        int len = *(int *)(v-1);
        free_dvector(v, 0, len);
    };
    py::capsule free_r(r, "free_r", dvlen);
    py::capsule free_z(z, "free_z", dvlen);

    return std::make_pair(
        py::array_t<double>(len, r, free_r), 
        py::array_t<double>(len, z, free_z)
    );
}

PYBIND11_MODULE(core, m) {
    // the order of definitions is important because
    // pybind11-stubgen can't handle forward declarations
    static Logger logger;                    // initialize the logfile
    m.doc() = "Python bindings for DipolEq"; // optional module docstring
    py::object __version__ = py::cast(VERSION_FULL);
    m.attr("__version__") = __version__;
    py::object __version_info__ = py::dict(
            "version"_a=VERSION,
            "version_full"_a=VERSION_FULL,
            "git_commit"_a=GIT_COMMIT,
            "git_short_hash"_a=GIT_SHORT_HASH,
            "git_distance"_a=GIT_DISTANCE
        );
    m.attr("__version_info__")=__version_info__;
    py::class_<TOKAMAK>(m, "Machine", "DipolEq machine object");   // forward declaration
    py::class_<DMatrixView>(m, "MatrixView", py::buffer_protocol(), "Matrix view class")
        .def_buffer(&DMatrixView::get_buffer_info)
    ;
    py::class_<DVectorView>(m, "VectorView", py::buffer_protocol(), "Vector view class")
        .def_buffer(&DVectorView::get_buffer_info)
    ;
    py::class_<IMatrixView>(m, "IMatrixView", py::buffer_protocol(), "Integer matrix view class")
        .def_buffer([](IMatrixView& self) -> py::buffer_info {
            return self.get_buffer_info();})
    ;

    py::enum_<ModelType>(m, "ModelType", py::arithmetic(), "Plasma model type")
        .value("Std", ModelType::Std)
        .value("IsoNoFlow", ModelType::IsoNoFlow)
        .value("IsoFlow", ModelType::IsoFlow)
        .value("AnisoNoFlow", ModelType::AnisoNoFlow)
        .value("AnisoFlow", ModelType::AnisoFlow)
        .value("DipoleStd", ModelType::DipoleStd)
        .value("DipoleIntStable", ModelType::DipoleIntStable)
        .value("DipoleStablePsiN", ModelType::DipoleStablePsiN)
    ;
    py::class_<CPlasmaModel>(m, "CPlasmaModel", "Plasma model base class")
        // FIXME
        // more work here... should make it possible to make 
        // our own python classes that inherit from CPlasmaModel
        // with overriden methods and generic getter and setters
        .def("update_model", &CPlasmaModel::UpdateModel, py::arg("Machine"), "Update the plasma model")
        .def("model_input", &CPlasmaModel::ModelInput, "Input model parameters")
    ;
    py::class_<PLASMA>(m, "Plasma", "plasma data and methods")
        .def(py::init(&new_Plasma), "Create Plasma")
        .def("init", &init_Plasma, "Initialize Plasma")
        .def("plasmaP", &PlasmaP, "Calculate plasma pressure")
        .def("plasmaPp", &PlasmaPp, "Calculate plasma Pprime")
        .def("plasmaG", &PlasmaG, "Calculate plasma G")
        .def("plasmaG2p", &PlasmaG2p, "Calculate plasma G2prime")
        .def_readwrite("Nsize", &PLASMA::Nsize, "Size of grid -1, should be factor of 2")
        .def_property("ModelType", [](PLASMA &self)
                      { return ModelType(self.ModelType); }, [](PLASMA &self, ModelType mtype)
                      { self.ModelType = (int)mtype; }, "Plasma model type, see ModelType enum")
        .def_property_readonly("Model", [](PLASMA &self) -> py::object
                               {
            if (self.Model == NULL) return py::none();
            return py::cast(self.Model); }, "The plasma model object")
        .def_readwrite("R0", &PLASMA::R0, "Reference major radius")
        .def_readwrite("Z0", &PLASMA::Z0, "Reference vertical position")
        .def_readwrite("B0", &PLASMA::B0, "Vacuum magnetic field at R0, Z0")
        .def_readwrite("Ip0", &PLASMA::Ip0, "initial plasma current")
        .def_readwrite("B0R0", &PLASMA::B0R0, ":math:`B_0 * R_0`")
        .def_readwrite("Jedge", &PLASMA::Jedge, "Edge current density")
        .def_readwrite("NumBndMomts", &PLASMA::NumBndMomts, "Number of moments of boundary to calculate")
        .def_readwrite("NumPsiPts", &PLASMA::NumPsiPts, "Number of normalized psi points to calculate")
        .def_readwrite("PsiXmax", &PLASMA::PsiXmax, "Outermost normalized Psi from 0.0 to 1.0 (say 0.995)")

        // old plasma terms
        .def_readwrite("G2pTerms", &PLASMA::G2pTerms, "Number of G^2' terms in polynomial plasma models")
        .def_readwrite("HTerms", &PLASMA::HTerms, "Number of H terms in polynomial plasma models")
        .def_readwrite("PpTerms", &PLASMA::PpTerms, "Number of p' terms in polynomial plasma models")
        .def_readwrite("RotTerms", &PLASMA::RotTerms, "Number of Rot terms in flow polynomial plasma models")
        .def_readwrite("SisoTerms", &PLASMA::SisoTerms, "Number of Siso terms in isotropic polynomial plasma models")
        .def_readwrite("SparTerms", &PLASMA::SparTerms, "Number of Spar terms in anisotropic polynomial plasma models")
        .def_readwrite("SperTerms", &PLASMA::SperTerms, "Number of Sper terms in anisotropic polynomial plasma models")
        .def_property_readonly("G2p", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.G2pTerms, self.G2p); }, "G^2' polynomial terms")
        .def_property_readonly("H", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.HTerms, self.H); }, "H polynomial terms")
        .def_property_readonly("Pp", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.PpTerms, self.Pp); }, "p' polynomial terms")
        .def_property_readonly("Rot", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.RotTerms, self.Rot); }, "Rotational flow terms")
        .def_property_readonly("Siso", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.SisoTerms, self.Siso); }, "Isotropic S polynomial terms")
        .def_property_readonly("Spar", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.SparTerms, self.Spar); }, "Anisotropic S parallel terms")
        .def_property_readonly("Sper", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.SperTerms, self.Sper); }, "Anisotropic S perpendicular terms")

        // results
        .def_readonly("Ip", &PLASMA::Ip, "Plasma current")
        .def_readonly("beta0", &PLASMA::beta0, "Vacuum toroidal beta at R0")
        .def_readonly("beta", &PLASMA::beta, "Average toroidal beta")
        .def_readonly("betap", &PLASMA::betap, "Poloidal beta")
        .def_readonly("li", &PLASMA::li, "Normalized internal inductance")
        .def_readonly("Ltotal", &PLASMA::Ltotal, "Total inductance")
        .def_readonly("mu", &PLASMA::mu, "Normalized diamagnetism")
        .def_readonly("Volume", &PLASMA::Volume, "Volume")
        .def_readonly("CrossSection", &PLASMA::CrossSection, "Cross section")
        .def_readonly("Perimeter", &PLASMA::Perimeter, "Perimeter")
        .def_readonly("Diamag", &PLASMA::Diamag, "Diamagnetism")
        .def_readonly("q0", &PLASMA::q0, "Central safety factor")
        .def_readonly("qCircular", &PLASMA::qCircular, "qCircular")
        .def_readonly("qStar", &PLASMA::qStar, "qStar")
        .def_readonly("RMagAxis", &PLASMA::XMagAxis, "Magnetic axis R")
        .def_readonly("ZMagAxis", &PLASMA::ZMagAxis, "Magnetic axis Z")
        .def_readonly("PsiMagAxis", &PLASMA::PsiMagAxis, "Psi at magnetic axis")
        .def_readonly("PsiAxis", &PLASMA::PsiAxis, "Psi at axis or inner separatrix")
        .def_readonly("PsiLim", &PLASMA::PsiLim, "Psi at plasma/vacuum boundary")
        .def_readonly("PsiFCFS", &PLASMA::PsiAxis, "Psi at axis or inner separatrix")
        .def_readonly("PsiLCFS", &PLASMA::PsiLim, "Psi at outer separatrix")
        .def_readonly("HalfWidth", &PLASMA::HalfWidth, "Half width")
        .def_readonly("Elongation", &PLASMA::Elongation, "Elongation")

        .def_readonly("ChiSqr", &PLASMA::ChiSqr, "Chi squared")
        .def_readonly("totKinEnergy", &PLASMA::TotKinEnergy, "Total kinetic energy")
        .def_readonly("totMagEnergy", &PLASMA::TotMagEnergy, "Total magnetic energy")
        .def_property_readonly("B2", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.B2); }, "B^2 on grid")
        .def_property_readonly("GradPsiR", [](PLASMA &self)
                               { return DMatrixView::create(self.Nsize, self.GradPsiX); }, "GradPsiR(self) -> MatrixView\n\n:math:`\\del\\psi/\\delR` on grid\n")
        .def_property_readonly("GradPsiZ", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.GradPsiZ); }, "∂Psi/∂Z on grid")
        .def_property_readonly("GradPsi2", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.GradPsi2); }, ":math:`|d\\psi/dR|^2 + |d\\psi/dZ|^2` on grid")
        .def_property_readonly("Bt", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Bt); }, "B_toroidal on grid")
        .def_property_readonly("G", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.G); }, ":math:`G = F/(B_0 R_0)` on grid")
        .def_property_readonly("Rho", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Rho); }, "Rho on grid")
        .def_property_readonly("Piso", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Piso); }, "Isotropic pressure on grid")
        .def_property_readonly("Ppar", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Ppar); }, "Anisotropic p_parallel pressure on grid")
        .def_property_readonly("Pper", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Pper); }, "Anistropic p_perpendicular pressure on grid")
        .def_property_readonly("Alpha", [](PLASMA &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Alpha); }, "Alpha = :math:`\\mu_0*(p_\\parallel - p_\\perp)/B^2` on grid")

        .def_property_readonly("PsiX_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.PsiX_pr); }, "Normalized Psi ordinate")
        .def_property_readonly("Psi_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Psi_pr); }, "Psi on normalized profile")
        .def_property_readonly("P_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.P_pr); }, "Pressure profile")
        .def_property_readonly("G_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.G_pr); }, "G = (F = B_t R)/(B_0 R_0) profile")
        .def_property_readonly("Pp_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Pp_pr); }, "p' profile")
        .def_property_readonly("G2p_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.G2p_pr); }, "G^2' profile")
        .def_property_readonly("q_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.q_pr); }, "q profile")
        .def_property_readonly("Volp_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Volp_pr); }, "dV/dPsi profile")
        .def_property_readonly("Vol_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Vol_pr); }, "V profile")
        .def_property_readonly("S_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.S_pr); }, "S = global shear of q profile")
        .def_property_readonly("B2_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.B2_pr); }, "<B^2> flux tube averaged profile")
        .def_property_readonly("Well_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Well_pr); }, "Magnetic well profile")
        .def_property_readonly("J_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.J_pr); }, "<J> flux tube averaged current density profile")
        .def_property_readonly("Beta_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.Beta_pr); }, "<β> flux tube averaged profile")
        .def_property_readonly("BetaMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.BetaMax_pr); }, "max(β) on flux tube profile")
        .def_property_readonly("RBetaMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.XBetaMax_pr); }, "R of max(β) profile")
        .def_property_readonly("ZBetaMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.ZBetaMax_pr); }, "Z of max(β) profile")
        .def_property_readonly("BBetaMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.BBetaMax_pr); }, "B of max(β) profile")
        .def_property_readonly("BMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.BMax_pr); }, "max(B) on flux tube profile")
        .def_property_readonly("RBMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.XBMax_pr); }, "R of max(B) profile")
        .def_property_readonly("ZBMax_pr", [](PLASMA &self) -> py::object
                               { return DVectorView::create(self.NumPsiPts, self.ZBMax_pr); }, "Z of max(B) profile");

    py::class_<PSIGRID>(m, "PsiGrid", "Magnetic flux grid data and methods")
        .def(py::init(&new_PsiGrid), "Create PsiGrid")
        .def("init", &init_PsiGrid, "Initialize PsiGrid")
        .def("go_PDE", &GoPDE, "Solve the PDE on the PsiGrid")
        .def("make_psi_symmetric", &MakePsiSymmetric, "Make Psi symmetric")
        .def("get_new_residual", &GetNewResidual, "Get new residual")
        .def("new_solution", &NewSolution, "Get new solution")
        .def("new_M_solution", &NewMSolution, "Get new M solution")
        .def("get_Psi", py::vectorize(&GetPsi), "Get Psi")
        .def("get_IsPlasma", &GetIsPlasma, "Get IsPlasma")
        .def("init_J", &InitJ, "Initialize J")
        .def_readwrite("Nsize", &PSIGRID::Nsize, "Size of grid -1, should be factor of 2")
        .def_readwrite("Symmetric", &PSIGRID::Symmetric, "Is the grid symmetric")
        .def_readwrite("MaxRes", &PSIGRID::MaxRes)
        .def_readwrite("PastMaxRes", &PSIGRID::PastMaxRes)
        .def_readwrite("Rmax", &PSIGRID::Xmax, "Maximum R of computational grid")
        .def_readwrite("Rmin", &PSIGRID::Xmin, "Minimum R of computational grid")
        .def_readwrite("Zmax", &PSIGRID::Zmax, "Maximum Z of computational grid")
        .def_readwrite("Zmin", &PSIGRID::Zmin, "Minimum Z of computational grid")
        .def_property_readonly("dr", [](PSIGRID &pg)
                               { return pg.dx; }, "R grid spacing")
        .def_property_readonly("dz", [](PSIGRID &pg)
                               { return pg.dz; }, "Z grid spacing")
        .def_readwrite("BoundError", &PSIGRID::BoundError)
        .def_readwrite("BoundThreshold", &PSIGRID::BoundThreshold)
        .def_readwrite("ResThreshold", &PSIGRID::ResThreshold)
        .def_readwrite("UnderRelax1", &PSIGRID::UnderRelax1)
        .def_readwrite("UnderRelax2", &PSIGRID::UnderRelax2)
        .def_readonly("PsiAxis", &PSIGRID::PsiAxis, "Psi at FCFS or Magnetic Axis")
        .def_readonly("PsiMagAxis", &PSIGRID::PsiMagAxis, "Psi at Magnetic Axis")
        .def_readonly("PsiLim", &PSIGRID::PsiLim, "Psi at plasma/vacuum boundary")
        .def_readonly("PsiFCFS", &PSIGRID::PsiAxis, "Psi at axis or inner separatrix")
        .def_readonly("PsiLCFS", &PSIGRID::PsiLim, "Psi at outer separatrix")
        .def_property_readonly("DelPsi", [](PSIGRID &pg)
                               { return pg.DelPsi; }, "Change in psi across plasma")
        .def_readwrite("RMagAxis", &PSIGRID::XMagAxis, "Magnetic axis R")
        .def_readwrite("ZMagAxis", &PSIGRID::ZMagAxis, "Magnetic axis Z")
        .def_property_readonly("R", [](PSIGRID &self) -> py::object
                               { return DVectorView::create(self.Nsize + 1, self.X); }, "R grid")
        .def_property_readonly("Z", [](PSIGRID &self) -> py::object
                               { return DVectorView::create(self.Nsize + 1, self.Z); }, "Z grid")
        .def_property_readonly("IsPlasma", [](PSIGRID &self) -> py::object
                               { return IMatrixView::create(self.Nsize, self.IsPlasma); }, "grid point is considered inside the plasma.")
        .def_property_readonly("Psi", [](PSIGRID &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Psi); }, "Magnetic flux on grid")
        .def_property_readonly("Current", [](PSIGRID &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Current); }, "Current density on grid * μ_0")
        .def_property_readonly("Residual", [](PSIGRID &self) -> py::object
                               { return DMatrixView::create(self.Nsize, self.Residual); }, "Current density residual on grid")
        .def("get_contour", &get_flux_contour, "Get a contour at normalized psi, returns r, z vectors");

    py::enum_<MeasType>(m, "MeasType", py::arithmetic(), "Measurement type")
        .value("unk", MeasType::unk)
        .value("bp", MeasType::bp)
        .value("press", MeasType::press)
        .value("pperp", MeasType::pperp)
        .value("ppar", MeasType::ppar)
        .value("flux", MeasType::flux)
        .value("saddle", MeasType::saddle)
        .value("circle", MeasType::circle)
        .value("coilcur", MeasType::coilcur)
        .value("plasmacur", MeasType::plasmacur)
        .value("bt", MeasType::bt)
        .value("diam", MeasType::diam)
        .value("bangle", MeasType::bangle)
        .value("flowt", MeasType::flowt)
        .value("flowp", MeasType::flowp)
        .value("ne", MeasType::ne)
        .value("Te", MeasType::Te)
        .value("Ti", MeasType::Ti)
        .value("Zeff", MeasType::Zeff)
        .value("rho", MeasType::rho)
        .value("rot", MeasType::rot)
        .value("ppsix", MeasType::ppsix)
        .value("bpangle", MeasType::bangle)
        .value("pnorm", MeasType::pnorm)
        .value("J0", MeasType::J0)
    ;

    py::enum_<CircleType>(m, "CircleType", py::arithmetic(), "Circle measurement type")
        .value("btcos", CircleType::btcos)
        .value("brsin", CircleType::brsin)
        .value("brcos", CircleType::brcos)
    ;

    py::class_<LIMITER>(m, "Limiter", "Limiter data")
        .def(py::init(&new_Limiter), "Create Limiter")
        .def_readwrite("R1", &LIMITER::X1)
        .def_readwrite("Z1", &LIMITER::Z1)
        .def_readwrite("R2", &LIMITER::X2)
        .def_readwrite("Z2", &LIMITER::Z2)
        .def_readwrite("Enabled", &LIMITER::Enabled, "1 = outer limimter, 0 = disabled, -1 = inner limiter")
        .def_property_readonly("PsiLim", [](LIMITER& l){return l.PsiMin;}, "Psi at limiting point")
        .def_property_readonly("RLim", [](LIMITER& l){return l.Xmin;}, "R at limiting point")
        .def_property_readonly("ZLim", [](LIMITER& l){return l.Zmin;}, "Z at limiting point")
        .def_property("Name", [](LIMITER& self) {return self.Name;},
            [](LIMITER& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(LIMITER::Name)-1);
            }, "Name of the limiter")
    ;

    py::class_<SUBCOIL>(m, "SubCoil", "Subcoil data")
        .def(py::init(&new_SubCoil), "Create SubCoil")
        .def_readwrite("R", &SUBCOIL::X)
        .def_readwrite("Z", &SUBCOIL::Z)
        .def_readwrite("Fraction", &SUBCOIL::CurrentFraction, "Fraction of current")
        .def_property("Name", [](SUBCOIL& self) {return self.Name;},
            [](SUBCOIL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(COIL::Name)-1);
            }, "Name of the subcoil")
    ;
    py::class_<ObjVecView<SUBCOIL>>(m, "SubCoils", "Subcoil vector view")
        .def("__getitem__", &ObjVecView<SUBCOIL>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<SUBCOIL>& self, size_t i, SUBCOIL* c) {
            if (i<0 || i>=self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<SUBCOIL>& self) {return self.size();})
        .def("new_subcoil", [](ObjVecView<SUBCOIL>& self) {
            SUBCOIL *c = new_SubCoil();
            return c;
        }, py::return_value_policy::reference_internal, "Add a new subcoil")
    ;

    py::class_<COIL>(m, "Coil", "Coil data and methods")
        .def(py::init(&new_Coil), "Create Coil")
        .def_readwrite("Enabled", &COIL::Enabled)
        .def_property("Name", [](COIL& self) {return self.Name;},
            [](COIL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(COIL::Name)-1);
            }, "Name of the coil")
        .def_property_readonly("SubCoils", [](COIL& self) -> py::object
            {return ObjVecView<SUBCOIL>::create(self.NumSubCoils, self.SubCoils);},
            "Return vector of SubCoils")
        .def_property("NumSubCoils", [](COIL& self) {return self.NumSubCoils;},
            &set_NumSubCoils, "Number of subcoils, setting will erase old subcoils")
        .def_readwrite("CoilCurrent", &COIL::CoilCurrent)
        .def_readwrite("R", &COIL::X, "Coil centroid R")
        .def_readwrite("dR", &COIL::dX, "Coil radial width, > 0 for automatic subcoil generation")
        .def_readwrite("Z", &COIL::Z, "Coil centroid Z")
        .def_readwrite("dZ", &COIL::dZ, "Coil vertical width")
        .def("compute_SubCoils", &compute_SubCoils, "Generate subcoils from coil parameters and grid")
    ;

    py::class_<SUBSHELL>(m, "SubShell", "Subshell data")
        .def(py::init(&new_SubShell), "Create SubShell")
        .def_property("Name", [](SUBSHELL& self) {return self.Name;},
            [](SUBSHELL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(SUBSHELL::Name)-1);
            }, "Name of the subshell")
        .def_readwrite("R", &SUBSHELL::X)
        .def_readwrite("Z", &SUBSHELL::Z)
        .def_readwrite("Radius", &SUBSHELL::Radius)
        .def_readwrite("Current", &SUBSHELL::Current)
    ;

    py::class_<ObjVecView<SUBSHELL>>(m, "SubShells", "Subshell vector view")
        .def("__getitem__", &ObjVecView<SUBSHELL>::operator[],
            py::return_value_policy::reference)
        .def("__len__", [](ObjVecView<SUBSHELL>& self) {return self.size();})
        // .def("__setitem__", &ObjVecView<SUBSHELL>::insert)
        // this doesn't work because the free operation can't be passed in
        // .def("new_subshell", [](ObjVecView<SUBSHELL>& self) { return new_SubShell();},
        //    py::return_value_policy::reference, "Add a new subshell")
    ;

    py::class_<SHELL>(m, "Shell", "Shell data and methods" )
        .def(py::init(&new_Shell), "Create Shell")
        .def_readwrite("Enabled", &SHELL::Enabled)
        .def_property("Name", [](SHELL& self) {return self.Name;},
            [](SHELL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(SHELL::Name)-1);
            }, "Name of the shell")
        .def_readonly("NumSubShells", &SHELL::NumSubShells)
        .def("set_NumSubShells", &set_SHELL_NumSubShells, "Set the number of subshells")
        .def_property_readonly("SubShells", [](SHELL& self) -> py::object
            {return ObjVecView<SUBSHELL>::create(self.NumSubShells, self.SubShells);},
             "Return vector of SubShells")
    ;

    py::class_<SEPARATRIX>(m, "Separatrix", "Separatrix data")
        .def(py::init(&new_Separatrix), "Create Separatrix")
        .def_property("Name", [](SEPARATRIX& self) {return self.Name;},
            [](SEPARATRIX& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(SEPARATRIX::Name)-1);
            }, "Name of the separatrix")
        .def_readwrite("Enabled", &SEPARATRIX::Enabled)
        .def_readwrite("IsSeparatrix", &SEPARATRIX::IsSeparatrix)
        .def_readwrite("R1", &SEPARATRIX::X1, "Box to search for separatrix")
        .def_readwrite("Z1", &SEPARATRIX::Z1, "Box to search for separatrix")
        .def_readwrite("R2", &SEPARATRIX::X2, "Box to search for separatrix")
        .def_readwrite("Z2", &SEPARATRIX::Z2, "Box to search for separatrix")
        .def_readwrite("RC", &SEPARATRIX::XC, "Center of plasma from Sep")
        .def_readwrite("ZC", &SEPARATRIX::ZC, "Center of plasma from Sep")
        .def_readonly("Psi", &SEPARATRIX::PsiSep, "Value of Psi at separatrix")
        .def_readonly("Rs", &SEPARATRIX::Xs, "R of separatrix")
        .def_readonly("Zs", &SEPARATRIX::Zs, "Z of separatrix")
    ;

    py::class_<MEAS>(m, "Measure", "Measurement data")
        .def(py::init(&new_Measure), "Create Measurement")
        .def_property("Name", [](MEAS& self) {return self.Name;},
            [](MEAS& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(MEAS::Name)-1);
            }, "Name of the measurement")
        .def_property("Type", [](MEAS& self){ return MeasType(self.mType);},
            [](MEAS& self, MeasType type){self.mType = (int) type;}, "Type of the measurement")
        .def_readwrite("R", &MEAS::X)
        .def_readwrite("Z", &MEAS::Z)
        .def_readwrite("Value", &MEAS::Value)
        .def_readwrite("StdDev", &MEAS::StdDev)
        .def_readwrite("Fit", &MEAS::Fit)
        .def_readwrite("Now", &MEAS::Now)
        // fixme: check for the right type
        .def_property("Angle", [](MEAS& self) {return self.parm.bp.Angle;},
            [](MEAS& self, double angle){self.parm.bp.Angle = angle;}, "Angle")
        .def_property("Radius", [](MEAS& self) {return self.parm.circle.Radius;},
            [](MEAS& self, double radius){self.parm.circle.Radius = radius;}, "Radius")
        .def_property("Number", [](MEAS& self) {return self.parm.circle.Number;},
            [](MEAS& self, int number){self.parm.circle.Number = number;}, "Number")
        .def_property("CircleType", [](MEAS& self) {return CircleType(self.parm.circle.CircleType);},
            [](MEAS& self, CircleType ct){self.parm.circle.CircleType = (int) ct;}, "CircleType")
        .def_property("CoilNum", [](MEAS& self) {return self.parm.coilcur.CoilNum;},
            [](MEAS& self, int num){self.parm.coilcur.CoilNum = num;}, "CoilNum")
        .def_property("R1", [](MEAS& self) {return self.parm.saddle.X1;},
            [](MEAS& self, double x1){self.parm.saddle.X1 = x1;}, "Saddle Loop R1")
        .def_property("Z1", [](MEAS& self) {return self.parm.saddle.Z1;},
            [](MEAS& self, double z1){self.parm.saddle.Z1 = z1;}, "Saddle Loop Z1")
        .def_property("R2", [](MEAS& self) {return self.parm.saddle.X2;},
            [](MEAS& self, double x2){self.parm.saddle.X2 = x2;}, "Saddle Loop R2")
        .def_property("Z2", [](MEAS& self) {return self.parm.saddle.Z2;},
            [](MEAS& self, double z2){self.parm.saddle.Z2 = z2;}, "Saddle Loop Z2")
    ;

    py::class_<ObjVecView<COIL>>(m, "Coils", "Coil vector view")
        .def("__getitem__", &ObjVecView<COIL>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", &ObjVecView<COIL>::insert)
        .def("__len__", [](ObjVecView<COIL>& self) {return self.size();})
        .def("new_Coil", [](ObjVecView<COIL>& self, int n) { return new_Coil(n);},
            py::return_value_policy::reference, "Add a new coil")
    ;
    py::class_<ObjVecView<LIMITER>>(m, "Limiters", "Limiter vector view")
        .def("__getitem__", &ObjVecView<LIMITER>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<LIMITER>& self, size_t i, LIMITER* c) {
            if (i<0 || i>=self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<LIMITER>& self) {return self.size();})
        .def("new_limiter", [](ObjVecView<LIMITER>& self) {
            LIMITER *c = new_Limiter();
            return c;
        }, py::return_value_policy::reference_internal, "Add a new limiter")
    ;
    py::class_<ObjVecView<SEPARATRIX>>(m, "Separatrices", "Separatrix vector view")
        .def("__getitem__", &ObjVecView<SEPARATRIX>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<SEPARATRIX>& self, size_t i, SEPARATRIX* c) {
            if (i<0 || i>=self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<SEPARATRIX>& self) {return self.size();})
        .def("new_separatrix", [](ObjVecView<SEPARATRIX>& self) {
            SEPARATRIX *c = new_Separatrix();
            return c;
        }, py::return_value_policy::reference, "Add a new separatrix")
    ;
    py::class_<ObjVecView<MEAS>>(m, "Measures", "Measurement vector view")
        .def("__getitem__", &ObjVecView<MEAS>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<MEAS>& self, size_t i, MEAS* c) {
            if (i<0 || i>=self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<MEAS>& self) {return self.size();})
        .def("new_meas", [](ObjVecView<MEAS>& self, MeasType mtype) {
            MEAS *c = new_Measure((int) mtype);
            return c;
        }, py::return_value_policy::reference, "Add a new measurement of type mtype")
    ;

    py::class_<ObjVecView<SHELL>>(m, "Shells", "Shell vector view")
        .def("__getitem__", &ObjVecView<SHELL>::operator[],
            py::return_value_policy::reference)
        .def("__setitem__", &ObjVecView<SHELL>::insert)
        .def("__len__", [](ObjVecView<SHELL>& self) {return self.size();})
        .def("new_shell", [](ObjVecView<SHELL>& self) { return new_Shell(0);},
            py::return_value_policy::reference, "Add a new shell")
    ;

    // and finally the Machine!
    py::class_<TOKAMAK>(m.attr("Machine"))
        .def(py::init(&new_Tokamak),  "Return new Machine struct",
            py::return_value_policy::reference)
        .def(py::init(&FileInput), "Initialize Machine with a file",
            py::return_value_policy::reference)
        .def("__del__", [](TOKAMAK& self) {free_Tokamak(&self, FALSE);},
            "Free the machine") // this causes double free errors
        .def("init", [](TOKAMAK& self) {init_Tokamak(&self);},
            "Initialize Machine")
        .def("set_coil_NumSubCoils",
                [](TOKAMAK& self, int i, int n) {
                    if (self.Coils[i] != NULL)
                        free_Coil(self.Coils[i], self.PsiGrid->Nsize);
                    self.Coils[i] = new_Coil(n);
                },"Set the number of subcoils for a coil")
        // This should be used, but when defined causes "double free" errors
        //.def("__del__", [](TOKAMAK& self) {printf("freeing %zu", (size_t)&self);
        //                                 free_Tokamak(&self);},
        //    "Free the machine")
        .def("set_start_time", &SetStartTime, "Set the start time")
        .def("set_stop_time", &SetStopTime, "Set the end time")
        .def("read_restart", [](TOKAMAK& self)  { ReadRestart(self.RSname, &self); }, "Read the restart file")
        .def("write_restart", [](TOKAMAK& self) { WriteRestart(self.RSname, &self); }, "Write the restart file")
        .def("find_shell_current", &Find_ShellCurrent, "Find shell current")
        .def("load_bndry_greens", &LoadBndryGreens, "Load boundary greens")
        .def("free_bndry_greens", &free_BndryGreens, "Load boundary greens")
        .def("psi_boundary", &PsiBoundary, "Calculate psi boundary")
        .def("add_coil_J", &AddCoilJ, "Add Coil currents")
        .def("add_shell_J", &AddShellJ, "Add Shell currents")
        .def("find_plasma_boundary", &PlasmaBoundary, "Find plasma boundary")
        .def("load_meas_greens", &LoadMeasGreens, "Load measurement greens")
        .def("free_meas_greens", &free_MeasGreens, "Free measurement greens")
        .def("least_squares", &LeastSquares, "Least squares")
        .def("get_plasma_parameters", &GetPlasmaParameters, "Get plasma parameters")
        .def("zero_J", &ZeroJ, "Zero J")
        .def("find_J", &FindJ, "Find J")
        .def_readwrite("MaxIterFixed", &TOKAMAK::MaxIterFixed)
        .def_readwrite("MaxIterFree", &TOKAMAK::MaxIterFree)
        .def_readwrite("IterFixed", &TOKAMAK::IterFixed)
        .def_readwrite("IterFree", &TOKAMAK::IterFree)
        .def_readwrite("LHGreenStatus", &TOKAMAK::LHGreenStatus)
        .def_readwrite("MGreenStatus", &TOKAMAK::MGreenStatus)
        .def_readwrite("SGreenStatus", &TOKAMAK::SGreenStatus)
        .def_readwrite("SInductStatus", &TOKAMAK::SInductStatus)
        .def_readwrite("RestartStatus", &TOKAMAK::RestartStatus)
        .def_readwrite("RestartUnkns", &TOKAMAK::RestartUnkns)
        .def_readwrite("VacuumOnly", &TOKAMAK::VacuumOnly)
        .def_readwrite("NumEqualEq", &TOKAMAK::NumEqualEq)
        .def_readwrite("NumMCarloEq", &TOKAMAK::NumMCarloEq)
        .def_readwrite("NumMCarloData", &TOKAMAK::NumMCarloData)
        .def_readwrite("MaxIterMCarlo", &TOKAMAK::MaxIterMCarlo)
        .def_readwrite("Confidence", &TOKAMAK::Confidence)
        .def_readonly("Start", &TOKAMAK::Start)
        .def_readonly("Stop",  &TOKAMAK::Stop)

        .def_readonly("PsiGrid", &TOKAMAK::PsiGrid,
            py::return_value_policy::reference_internal)
        .def_readonly("Plasma", &TOKAMAK::Plasma,
            py::return_value_policy::reference_internal)

        .def_property("Name",[](TOKAMAK& self) {return self.Info;},
            [](TOKAMAK& self, std::string name){
            strncpy(self.Name, name.c_str(), sizeof(TOKAMAK::Name)-1);
            }, "Set the name")
        .def_property("Info", [](TOKAMAK& self) {return self.Info;},
            [](TOKAMAK& self, std::string info){
                strncpy(self.Info, info.c_str(), sizeof(TOKAMAK::Info)-1);
            }, "Info about the machine")
        .def_property("Iname", [](TOKAMAK& self) {return self.Iname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.Iname, name.c_str(), sizeof(TOKAMAK::Iname)-1);
            }, "Input name")
        .def_property("Oname", [](TOKAMAK& self) {return self.Oname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.Oname, name.c_str(), sizeof(TOKAMAK::Oname)-1);
            }, "Output name")
        .def_property("Iname", [](TOKAMAK& self) {return self.Iname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.Iname, name.c_str(), sizeof(TOKAMAK::Iname)-1);
            }, "Input name")
        .def_property("LHname", [](TOKAMAK& self) {return self.LHname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.LHname, name.c_str(), sizeof(TOKAMAK::LHname)-1);
            }, "LH name")
        .def_property("MGname", [](TOKAMAK& self) {return self.MGname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.MGname, name.c_str(), sizeof(TOKAMAK::MGname)-1);
            }, "MG name")
        .def_property("SGname", [](TOKAMAK& self) {return self.SGname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.SGname, name.c_str(), sizeof(TOKAMAK::SGname)-1);
            }, "SG name")
        .def_property("SMname", [](TOKAMAK& self) {return self.SMname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.SMname, name.c_str(), sizeof(TOKAMAK::SMname)-1);
            }, "SM name")
        .def_property("RSname", [](TOKAMAK& self) {return self.RSname;},
            [](TOKAMAK& self, std::string name){
                strncpy(self.RSname, name.c_str(), sizeof(TOKAMAK::RSname)-1);
            }, "RS name")
        .def_property_readonly("Coils", [](TOKAMAK& self) -> py::object 
            {return ObjVecView<COIL>::create(self.NumCoils, self.Coils, &self, free_COIL);},
            "Array of COILS")
        .def_property("NumCoils", [](TOKAMAK& self) {return self.NumCoils;},
                    &set_NumCoils, "Number of coils, setting will erase old coils")
        .def("set_NumSubcoils", [](TOKAMAK& self, int i, int n) {
                if (self.Coils[i] != NULL)
                    free_Coil(self.Coils[i], self.PsiGrid->Nsize);
                self.Coils[i] = new_Coil(n);
            }, "Set the number of subcoils for a coil")
        .def_property_readonly("Limiters", [](TOKAMAK& self) -> py::object 
            { return ObjVecView<LIMITER>::create(self.NumLimiters, self.Limiters);},
            "Get the limiters")
        .def_property("NumLimiters", [](TOKAMAK& self) {return self.NumLimiters;},
            &set_NumLimiters, "Number of limiters, setting will erase old limiters")
        .def_property_readonly("Shells", [](TOKAMAK& self) -> py::object {
            return ObjVecView<SHELL>::create(self.NumShells, self.Shells, &self, free_SHELL);
            }, "Get the shells")
        .def_property("NumShells", [](TOKAMAK& self) {return self.NumShells;},
            &set_NumShells, "Number of shells, setting will erase old shells")
        .def("set_NumSubShells", &set_TOKAMAK_NumSubShells,
            "Set the number of subshells of shell i to n")
        .def_property_readonly("Measures", [](TOKAMAK& self) -> py::object {
            return ObjVecView<MEAS>::create(self.NumMeasures, self.Measures, &self, free_MEAS);
            }, "Get the measurements")
        .def_property("NumMeasures", [](TOKAMAK& self) {return self.NumMeasures;},
            &set_NumMeasures, "Number of measurements, setting will erase old measurements")
        .def_property_readonly("Seps", [](TOKAMAK& self) -> py::object {
            return ObjVecView<SEPARATRIX>::create(self.NumSeps, self.Seps);
            }, "Get the separatrixes")
        .def_property("NumSeps", [](TOKAMAK& self) {return self.NumSeps;},
            &set_NumSeps, "Number of separatrixes, setting will erase old separatrixes")
    ;

}
