#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
#include "PsiBoundary.h"
#include "PlasmaBoundary.h"
#include "Find_ShellCurrent.h"
#include "FindMeasFit.h"


namespace py = pybind11;

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

class DMatrixView {
public:
    DMatrixView(size_t nsize, double ** data) : m_size(nsize), m_data(data) {};
    double& operator()(size_t i, size_t j) { return m_data[i][j]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data[0],                             /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            2,                                      /* Number of dimensions */
            { m_size, m_size },                 /* Buffer dimensions */
            { sizeof(double) * m_size,            /* Strides (in bytes) for each index */
              sizeof(double) }
        );
    }
private:
    size_t m_size;
    double ** m_data;
};

class IMatrixView {
public:
    IMatrixView(size_t nsize, int ** data) : m_size(nsize), m_data(data) {};
    int& operator()(size_t i, size_t j) { return m_data[i][j]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data[0],                             /* Pointer to buffer */
            sizeof(int),                          /* Size of one scalar */
            py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
            2,                                      /* Number of dimensions */
            { m_size, m_size },                 /* Buffer dimensions */
            { sizeof(int) * m_size,            /* Strides (in bytes) for each index */
              sizeof(int) }
        );
    }
private:
    size_t m_size;
    int ** m_data;
};

class DVectorView {
public:
    DVectorView(size_t nsize, double * data) : m_size(nsize), m_data(data) {};
    double& operator[](size_t i) { return m_data[i]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data,                             /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            1,                                      /* Number of dimensions */
            { m_size },                 /* Buffer dimensions */
            { sizeof(double) }
        );
    }
private:
    size_t m_size;
    double * m_data;
};

extern "C" {
    extern FILE *LogFile;
};

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

enum class ModelType {
    Std = Plasma_Std,
    IsoNoFlow = Plasma_IsoNoFlow,
    IsoFlow = Plasma_IsoFlow,
    AnisoNoFlow = Plasma_AnisoNoFlow,
    AnisoFlow = Plasma_AnisoFlow,
    DipoleStd = Plasma_DipoleStd,
    DipoleIntStable = Plasma_DipoleIntStable
};

PYBIND11_MODULE(_dipoleq, m) {
    static Logger logger;                    // initialize the logfile
    m.doc() = "Python bindings for DipolEq"; // optional module docstring
    
    py::class_<DMatrixView>(m, "MatrixView", py::buffer_protocol())
        .def_buffer(&DMatrixView::get_buffer_info)
    ;
    py::class_<DVectorView>(m, "VectorView", py::buffer_protocol())
        .def_buffer(&DVectorView::get_buffer_info)
    ;
    py::class_<IMatrixView>(m, "IMatrixView", py::buffer_protocol())
        .def_buffer(&IMatrixView::get_buffer_info)
    ;
    py::class_<TOKAMAK>(m, "TOKAMAK")
        .def(py::init(&new_Tokamak), "Return new TOKAMAK struct")
        .def("init", [](TOKAMAK& self) {init_Tokamak(&self);},
            "Initialize TOKAMAK")
        .def(py::init(&FileInput), "Initialize TOKAMAK with a file")
        .def("__del__", &free_Tokamak, "Free TOKAMAK and what is below it")
        .def("set_start_time", &SetStartTime, "Set the start time")
        .def("set_stop_time", &SetStopTime, "Set the end time")
        .def("find_shell_current", &Find_ShellCurrent, "Find shell current")
        .def("load_bndry_greens", &LoadBndryGreens, "Load boundary greens")
        .def("psi_boundary", &PsiBoundary, "Calculate psi boundary")
        .def("add_coil_j", &AddCoilJ, "Add Coil currents")
        .def("add_shell_j", &AddShellJ, "Add Shell currents")
        .def("find_plasma_boundary", &PlasmaBoundary, "Find plasma boundary")
        .def("load_meas_greens", &LoadMeasGreens, "Load measurement greens")
        .def("least_squares", &LeastSquares, "Least squares")
        .def("free_meas_greens", &free_MeasGreens, "Free measurement greens")
        .def("find_shell_current", &Find_ShellCurrent, "Find shell current")
        .def("zero_j", &ZeroJ, "Zero J")
        .def("find_j", &FindJ, "Find J")
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

        .def_readwrite("NumCoils", &TOKAMAK::NumCoils)
        .def_readwrite("NumShells", &TOKAMAK::NumShells)
        .def_readwrite("NumLimiters", &TOKAMAK::NumLimiters)
        .def_readwrite("NumSeps", &TOKAMAK::NumSeps)
        .def_readwrite("NumMeasures", &TOKAMAK::NumMeasures)
        .def_readwrite("NumUnkns", &TOKAMAK::NumUnkns)        

        .def_readwrite("Confidence", &TOKAMAK::Confidence)
        .def_readonly("PsiGrid", &TOKAMAK::PsiGrid, py::return_value_policy::reference)

        .def_readonly("Name", &TOKAMAK::Name)
        .def_readonly("Info", &TOKAMAK::Info)
        .def_readonly("Iname", &TOKAMAK::Iname)
        .def_readonly("Oname", &TOKAMAK::Oname)
        .def_readonly("LHname", &TOKAMAK::LHname)
        .def_readonly("MGname", &TOKAMAK::MGname)
        .def_readonly("SGname", &TOKAMAK::SGname)
        .def_readonly("SMname", &TOKAMAK::SMname)
        .def_readonly("RSname", &TOKAMAK::RSname)
        .def_readonly("Start", &TOKAMAK::Start)
        .def_readonly("Stop", &TOKAMAK::Stop)
    ;

    py::class_<PSIGRID>(m, "PSIGRID")
        .def(py::init(&new_PsiGrid), "Create PSIGRID")
        .def("init", &init_PsiGrid, "Initialize PSIGRID")
        .def("go_PDE", &GoPDE, "Solve the PDE on the PSIGRID")
        .def("makePsiSymmetric", &MakePsiSymmetric, "Make Psi symmetric")
        .def("getNewResidual", &GetNewResidual, "Get new residual")
        .def("newSolution", &NewSolution, "Get new solution")
        .def("newMSolution", &NewMSolution, "Get new M solution")
        .def("getPsi", &GetPsi, "Get Psi")
        .def("getIsPlasma", &GetIsPlasma, "Get IsPlasma")
        .def_readwrite("Nsize", &PSIGRID::Nsize)
        .def_readwrite("Symmetric", &PSIGRID::Symmetric)
        .def_readwrite("MaxRes", &PSIGRID::MaxRes)
        .def_readwrite("PastMaxRes", &PSIGRID::PastMaxRes)
        .def_readwrite("Xmax", &PSIGRID::Xmax)
        .def_readwrite("Xmin", &PSIGRID::Xmin)
        .def_readwrite("Zmax", &PSIGRID::Zmax)
        .def_readwrite("Zmin", &PSIGRID::Zmin)
        .def_readwrite("dx", &PSIGRID::dx)
        .def_readwrite("dz", &PSIGRID::dz)
        .def_readwrite("BoundError", &PSIGRID::BoundError)
        .def_readwrite("BoundThreshold", &PSIGRID::BoundThreshold)
        .def_readwrite("ResThreshold", &PSIGRID::ResThreshold)
        .def_readwrite("UnderRelax1", &PSIGRID::UnderRelax1)
        .def_readwrite("UnderRelax2", &PSIGRID::UnderRelax2)
        .def_readwrite("PsiAxis", &PSIGRID::PsiAxis, "Psi at FCFS or Magnetic Axis")
        .def_readwrite("PsiMagAxis", &PSIGRID::PsiMagAxis, "Psi at Magnetic Axis")
        .def_readwrite("PsiLim", &PSIGRID::PsiLim, "Psi at plasma/vacuum boundary")
        .def_readwrite("DelPsi", &PSIGRID::DelPsi)
        .def_readwrite("XMagAxis", &PSIGRID::XMagAxis)
        .def_readwrite("ZMagAxis", &PSIGRID::ZMagAxis)
        .def("X", [](PSIGRID& self) {return DVectorView(self.Nsize, self.X);})
        .def("Z", [](PSIGRID& self) {return DVectorView(self.Nsize, self.Z);})
        .def("IsPlasma", [](PSIGRID& self) {return IMatrixView(self.Nsize, self.IsPlasma);})
        .def("Psi", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Psi);})
        .def("Current", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Current);})
        .def("Residual", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Residual);})
    ;    

    py::class_<CPlasmaModel>(m, "CPlasmaModel")
        .def("updateModel", &CPlasmaModel::UpdateModel, "Update the plasma model")
    ;

    py::class_<PLASMA>(m, "PLASMA")
        .def(py::init(&new_Plasma), "Create PLASMA")
        .def("init", &init_Plasma, "Initialize PLASMA")
        .def("plasmaP", &PlasmaP, "Calculate plasma pressure")
        .def("plasmaPp", &PlasmaPp, "Calculate plasma Pprime")
        .def("plasmaG", &PlasmaG, "Calculate plasma G")
        .def("plasmaG2p", &PlasmaG2p, "Calculate plasma G2prime")
        .def_readwrite("Nsize", &PLASMA::Nsize)
        .def_readwrite("ModelType", &PLASMA::ModelType)
        .def("Model", [](PLASMA& self) {return self.Model;}, py::return_value_policy::reference)
        .def_readwrite("R0", &PLASMA::R0, "Reference major radius")
        .def_readwrite("Z0", &PLASMA::Z0, "Reference vertical position")
        .def_readwrite("B0", &PLASMA::B0, "Vacuum magnetic field at R0, Z0")
        .def_readwrite("Ip0", &PLASMA::Ip0, "initial plasma current")
        .def_readwrite("B0R0", &PLASMA::B0R0, "B0 * R0")
        .def_readwrite("Jedge", &PLASMA::Jedge, "Edge current density")

        .def_readwrite("NumBndMomts", &PLASMA::NumBndMomts)
        .def_readwrite("NumPsiPts", &PLASMA::NumPsiPts)
        .def_readwrite("PsiXmax", &PLASMA::PsiXmax, "Outermost normalized Psi from 0.0 to 1.0")

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
        .def_readonly("XMagAxis", &PLASMA::XMagAxis, "Magnetic axis X")
        .def_readonly("ZMagAxis", &PLASMA::ZMagAxis, "Magnetic axis Z")
        .def_readonly("PsiMagAxis", &PLASMA::PsiMagAxis, "Psi at magnetic axis")
        .def_readonly("PsiAxis", &PLASMA::PsiAxis, "Psi at axis")
        .def_readonly("PsiLim", &PLASMA::PsiLim, "Psi at plasma/vacuum boundary")
        .def_readonly("HalfWidth", &PLASMA::HalfWidth, "Half width")
        .def_readonly("Elongation", &PLASMA::Elongation, "Elongation")

        .def_readonly("ChiSqr", &PLASMA::ChiSqr, "Chi squared")
        .def_readonly("totKinEnergy", &PLASMA::TotKinEnergy, "Total kinetic energy")
        .def_readonly("totMagEnergy", &PLASMA::TotMagEnergy, "Total agnetic energy")
        .def("B2", [](PLASMA& self) {return DMatrixView(self.Nsize, self.B2);})
        .def("GradPsiX", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsiX);})
        .def("GradPsiZ", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsiZ);})
        .def("GradPsi2", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsi2);})
        .def("Bt", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Bt);})
        .def("G", [](PLASMA& self) {return DMatrixView(self.Nsize, self.G);})
        .def("Rho", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Rho);})
        .def("Piso", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Piso);})
        .def("Ppar", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Ppar);})
        .def("Pper", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Pper);})
        .def("Alpha", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Alpha);})
    ;

    py::enum_<ModelType>(m, "ModelType")
        .value("Std", ModelType::Std)
        .value("IsoNoFlow", ModelType::IsoNoFlow)   
        .value("IsoFlow", ModelType::IsoFlow)
        .value("AnisoNoFlow", ModelType::AnisoNoFlow)
        .value("AnisoFlow", ModelType::AnisoFlow)
        .value("DipoleStd", ModelType::DipoleStd)
        .value("DipoleIntStable", ModelType::DipoleIntStable)
        .export_values()
    ;

    py::class_<SUBSHELL>(m, "SUBSHELL")
        .def(py::init(&new_SubShell), "Create SUBSHELL")
        .def_readonly("Name", &SUBSHELL::Name)
        .def_readwrite("X", &SUBSHELL::X)
        .def_readwrite("Z", &SUBSHELL::Z)
        .def_readwrite("Radius", &SUBSHELL::Radius)
        .def_readwrite("Current", &SUBSHELL::Current)
    ;

    py::class_<SHELL>(m, "SHELL")
        .def(py::init(&new_Shell), "Create SHELL")
        .def("add_subshell", &add_SubShell, "Add a subshell")
        .def_readwrite("Enabled", &SHELL::Enabled)
        .def_readonly("Name", &SHELL::Name)
        .def_readwrite("NumSubShells", &SHELL::NumSubShells)
        .def("SubShells", [](SHELL& self, int i) {return self.SubShells[i];}, py::return_value_policy::reference)
    ;
}
