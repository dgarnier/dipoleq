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
#include "Restart.h"


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
    DMatrixView(size_t nsize, double ** data) : m_size(nsize+1), m_data(data) {};
    double& operator()(size_t i, size_t j) { return m_data[i][j]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data[0],                             /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            2,                                      /* Number of dimensions */
            { m_size, m_size },                 /* Buffer dimensions */
            { sizeof(double) ,            /* Strides (in bytes) for each index */
              sizeof(double) * m_size }   /* nrutil is fortran order */
        );
    }
private:
    size_t m_size;
    double ** m_data;
};

class IMatrixView {
public:
    IMatrixView(size_t nsize, int ** data) : m_size(nsize+1), m_data(data) {};
    int& operator()(size_t i, size_t j) { return m_data[i][j]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data[0],                             /* Pointer to buffer */
            sizeof(int),                          /* Size of one scalar */
            py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
            2,                                      /* Number of dimensions */
            { m_size, m_size },                 /* Buffer dimensions */
            { sizeof(int) ,            /* Strides (in bytes) for each index */
              sizeof(int) * m_size}     /* nrutil is fortran order */
        );
    }
private:
    size_t m_size;
    int ** m_data;
};

class DVectorView {
public:
    DVectorView(size_t nsize, double * data) : m_size(nsize+1), m_data(data) {};
    double& operator[](size_t i) { return m_data[i]; };
    py::buffer_info get_buffer_info() {
        return py::buffer_info(
            m_data,                                  /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            1,                                      /* Number of dimensions */
            { m_size },                             /* Buffer dimensions */
            { sizeof(double) }
        );
    }
private:
    size_t m_size;
    double * m_data;
};

template <typename T> class ObjVecView  
{
    public:
    ObjVecView(size_t nsize, T ** data) : m_size(nsize), m_data(data) {};
    ObjVecView(size_t nsize, T ** data, TOKAMAK * mach, void (*objfree)(T *, TOKAMAK *)) 
        : m_size(nsize), m_data(data), m_machine(mach), m_free(objfree) {};

    T*& operator[](size_t i) { 
        if (i<0 || i>m_size)
            throw std::out_of_range("Index out of range");
        return m_data[i]; };
    size_t size() {return m_size;};

    void insert(size_t i, T* obj) {
        if (i<0 || i>m_size)
            throw std::out_of_range("Index out of range");
        if (m_data[i] != NULL) {
            if (m_free != NULL) {
                m_free(m_data[i], m_machine);
            } else {
                free(m_data[i]);
            }
        }
        m_data[i] = obj;
    }
    
    private:
    TOKAMAK *m_machine = NULL;
    void (*m_free)(T *, TOKAMAK *) = NULL;
    size_t m_size;
    T ** m_data;
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


PYBIND11_MODULE(_c, m) {
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
    py::class_<ObjVecView<COIL>>(m, "COILS")
        .def("__getitem__", &ObjVecView<COIL>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", &ObjVecView<COIL>::insert)
        .def("__len__", [](ObjVecView<COIL>& self) {return self.size();})
        .def("new_Coil", [](ObjVecView<COIL>& self, int n) { return new_Coil(n);}, 
            py::return_value_policy::reference, "Add a new coil")
    ;
    py::class_<ObjVecView<SUBCOIL>>(m, "SUBCOILS")
        .def("__getitem__", &ObjVecView<SUBCOIL>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<SUBCOIL>& self, size_t i, SUBCOIL* c) {
            if (i<0 || i>self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<SUBCOIL>& self) {return self.size();})
        .def("new_subcoil", [](ObjVecView<SUBCOIL>& self) {
            SUBCOIL *c = new_SubCoil();
            return c;
        }, py::return_value_policy::reference_internal, "Add a new subcoil")
    ;
    py::class_<ObjVecView<LIMITER>>(m, "LIMITERS")
        .def("__getitem__", &ObjVecView<LIMITER>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<LIMITER>& self, size_t i, LIMITER* c) {
            if (i<0 || i>self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<LIMITER>& self) {return self.size();})
        .def("new_limiter", [](ObjVecView<LIMITER>& self) {
            LIMITER *c = new_Limiter();
            return c;
        }, py::return_value_policy::reference_internal, "Add a new limiter")    
    ;
    py::class_<ObjVecView<SEPARATRIX>>(m, "SEPARATRICIES")
        .def("__getitem__", &ObjVecView<SEPARATRIX>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<SEPARATRIX>& self, size_t i, SEPARATRIX* c) {
            if (i<0 || i>self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<SEPARATRIX>& self) {return self.size();})
        .def("new_separatrix", [](ObjVecView<SEPARATRIX>& self) {
            SEPARATRIX *c = new_Separatrix();
            return c;
        }, py::return_value_policy::reference, "Add a new separatrix")
    ;
    py::class_<ObjVecView<MEAS>>(m, "MEASUREMENTS")
        .def("__getitem__", &ObjVecView<MEAS>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", [](ObjVecView<MEAS>& self, size_t i, MEAS* c) {
            if (i<0 || i>self.size())
                throw std::out_of_range("Index out of range");
            self[i] = c;
        })
        .def("__len__", [](ObjVecView<MEAS>& self) {return self.size();})
        .def("new_meas", [](ObjVecView<MEAS>& self, int mtype) {
            MEAS *c = new_Measure(0);
            return c;
        }, py::return_value_policy::reference, "Add a new measurement of type mtype")
    ;

    py::class_<ObjVecView<SHELL>>(m, "SHELLS")
        .def("__getitem__", &ObjVecView<SHELL>::operator[], 
            py::return_value_policy::reference)
        .def("__setitem__", &ObjVecView<SHELL>::insert)
        .def("__len__", [](ObjVecView<SHELL>& self) {return self.size();})
        .def("new_shell", [](ObjVecView<SHELL>& self) { return new_Shell(0);}, 
            py::return_value_policy::reference, "Add a new shell")
    ;

    py::class_<ObjVecView<SUBSHELL>>(m, "SUBSHELLS")
        .def("__getitem__", &ObjVecView<SUBSHELL>::operator[], 
            py::return_value_policy::reference)
        .def("__len__", [](ObjVecView<SUBSHELL>& self) {return self.size();})
        // .def("__setitem__", &ObjVecView<SUBSHELL>::insert)
        // this doesn't work because the free operation can't be passed in
        // .def("new_subshell", [](ObjVecView<SUBSHELL>& self) { return new_SubShell();}, 
        //    py::return_value_policy::reference, "Add a new subshell")
    ;

    py::class_<TOKAMAK>(m, "MACHINE")
        .def(py::init(&new_Tokamak), "Return new MACHINE struct")
        .def("init", [](TOKAMAK& self) {init_Tokamak(&self);},
            "Initialize MACHINE")
        .def("set_coil_NumSubCoils", 
                [](TOKAMAK& self, int i, int n) {
                    if (self.Coils[i] != NULL)
                        free_Coil(self.Coils[i], self.PsiGrid->Nsize);
                    self.Coils[i] = new_Coil(n);
                },"Set the number of subcoils for a coil")
        .def(py::init(&FileInput), "Initialize MACHINE with a file")
        .def("__del__", &free_Tokamak, "Free MACHINE and what is below it")
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
        .def("find_shell_current", &Find_ShellCurrent, "Find shell current")
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
        .def_property_readonly("Coils", [](TOKAMAK& self) {
            return ObjVecView<COIL>(self.NumCoils, self.Coils, &self, free_COIL);}, 
            "Array of COILS")
        .def_property("NumCoils", [](TOKAMAK& self) {return self.NumCoils;}, 
                    &set_NumCoils, "Number of coils, setting will erase old coils")
        .def("set_NumSubcoils", [](TOKAMAK& self, int i, int n) {
                if (self.Coils[i] != NULL)
                    free_Coil(self.Coils[i], self.PsiGrid->Nsize);
                self.Coils[i] = new_Coil(n);
            }, "Set the number of subcoils for a coil")
        .def_property_readonly("Limiters", [](TOKAMAK& self) {
            return ObjVecView<LIMITER>(self.NumLimiters, self.Limiters);
            }, "Get the limiters")
        .def_property("NumLimiters", [](TOKAMAK& self) {return self.NumLimiters;}, 
            &set_NumLimiters, "Number of limiters, setting will erase old limiters")
        .def_property_readonly("Shells", [](TOKAMAK& self) {
            return ObjVecView<SHELL>(self.NumShells, self.Shells, &self, free_SHELL);
            }, "Get the shells")
        .def_property("NumShells", [](TOKAMAK& self) {return self.NumShells;}, 
            &set_NumShells, "Number of shells, setting will erase old shells")
        .def("set_NumSubShells", &set_TOKAMAK_NumSubShells, 
            "Set the number of subshells of shell i to n")
        .def_property_readonly("Measures", [](TOKAMAK& self) {
            return ObjVecView<MEAS>(self.NumMeasures, self.Measures, &self, free_MEAS);
            }, "Get the measurements")
        .def_property("NumMeasures", [](TOKAMAK& self) {return self.NumMeasures;}, 
            &set_NumMeasures, "Number of measurements, setting will erase old measurements")
        .def_property_readonly("Seps", [](TOKAMAK& self) {
            return ObjVecView<SEPARATRIX>(self.NumSeps, self.Seps);
            }, "Get the separatrixes")
        .def_property("NumSeps", [](TOKAMAK& self) {return self.NumSeps;},
            &set_NumSeps, "Number of separatrixes, setting will erase old separatrixes")
    ;

    py::class_<PSIGRID>(m, "PSIGRID")
        .def(py::init(&new_PsiGrid), "Create PSIGRID")
        .def("init", &init_PsiGrid, "Initialize PSIGRID")
        .def("go_PDE", &GoPDE, "Solve the PDE on the PSIGRID")
        .def("make_psi_symmetric", &MakePsiSymmetric, "Make Psi symmetric")
        .def("get_new_residual", &GetNewResidual, "Get new residual")
        .def("new_solution", &NewSolution, "Get new solution")
        .def("new_M_solution", &NewMSolution, "Get new M solution")
        .def("get_Psi", py::vectorize(&GetPsi), "Get Psi")
        .def("get_IsPlasma", &GetIsPlasma, "Get IsPlasma")
        .def("init_J", &InitJ, "Initialize J")
        .def_readwrite("Nsize", &PSIGRID::Nsize)
        .def_readwrite("Symmetric", &PSIGRID::Symmetric)
        .def_readwrite("MaxRes", &PSIGRID::MaxRes)
        .def_readwrite("PastMaxRes", &PSIGRID::PastMaxRes)
        .def_readwrite("Rmax", &PSIGRID::Xmax)
        .def_readwrite("Rmin", &PSIGRID::Xmin)
        .def_readwrite("Zmax", &PSIGRID::Zmax)
        .def_readwrite("Zmin", &PSIGRID::Zmin)
        .def_readwrite("dr", &PSIGRID::dx)
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
        .def_readwrite("RMagAxis", &PSIGRID::XMagAxis)
        .def_readwrite("ZMagAxis", &PSIGRID::ZMagAxis)
        .def_property_readonly("R", [](PSIGRID& self) {return DVectorView(self.Nsize, self.X);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Z", [](PSIGRID& self) {return DVectorView(self.Nsize, self.Z);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("IsPlasma", [](PSIGRID& self) {return IMatrixView(self.Nsize, self.IsPlasma);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Psi", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Psi);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Current", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Current);}, 
            py::return_value_policy::reference_internal)
        .def_property_readonly("Residual", [](PSIGRID& self) {return DMatrixView(self.Nsize, self.Residual);},
            py::return_value_policy::reference_internal)
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
        .def_readonly("RMagAxis", &PLASMA::XMagAxis, "Magnetic axis R")
        .def_readonly("ZMagAxis", &PLASMA::ZMagAxis, "Magnetic axis Z")
        .def_readonly("PsiMagAxis", &PLASMA::PsiMagAxis, "Psi at magnetic axis")
        .def_readonly("PsiAxis", &PLASMA::PsiAxis, "Psi at axis")
        .def_readonly("PsiLim", &PLASMA::PsiLim, "Psi at plasma/vacuum boundary")
        .def_readonly("HalfWidth", &PLASMA::HalfWidth, "Half width")
        .def_readonly("Elongation", &PLASMA::Elongation, "Elongation")

        .def_readonly("ChiSqr", &PLASMA::ChiSqr, "Chi squared")
        .def_readonly("totKinEnergy", &PLASMA::TotKinEnergy, "Total kinetic energy")
        .def_readonly("totMagEnergy", &PLASMA::TotMagEnergy, "Total agnetic energy")
        .def_property_readonly("B2", [](PLASMA& self) {return DMatrixView(self.Nsize, self.B2);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("GradPsiX", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsiX);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("GradPsiZ", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsiZ);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("GradPsi2", [](PLASMA& self) {return DMatrixView(self.Nsize, self.GradPsi2);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Bt", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Bt);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("G", [](PLASMA& self) {return DMatrixView(self.Nsize, self.G);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Rho", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Rho);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Piso", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Piso);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Ppar", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Ppar);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Pper", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Pper);},
            py::return_value_policy::reference_internal)
        .def_property_readonly("Alpha", [](PLASMA& self) {return DMatrixView(self.Nsize, self.Alpha);},
            py::return_value_policy::reference_internal)
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

    py::class_<LIMITER>(m, "LIMITER")
        .def(py::init(&new_Limiter), "Create LIMITER")
        .def_readwrite("R1", &LIMITER::X1)
        .def_readwrite("Z1", &LIMITER::Z1)
        .def_readwrite("R2", &LIMITER::X2)
        .def_readwrite("Z2", &LIMITER::Z2)
        .def_readwrite("Enabled", &LIMITER::Enabled)
        .def_readwrite("PsiMin", &LIMITER::PsiMin)
        .def_readwrite("Rmin", &LIMITER::Xmin)
        .def_readwrite("Zmin", &LIMITER::Zmin)
        .def_property("Name", [](LIMITER& self) {return self.Name;},
            [](LIMITER& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(LIMITER::Name)-1);
            }, "Name of the limiter")
    ;

    py::class_<SUBCOIL>(m, "SUBCOIL")
        .def(py::init(&new_SubCoil), "Create SUBCOIL")
        .def_readwrite("R", &SUBCOIL::X)
        .def_readwrite("Z", &SUBCOIL::Z)
        .def_readwrite("Fraction", &SUBCOIL::CurrentFraction, "Fraction of current")
        .def_property("Name", [](SUBCOIL& self) {return self.Name;},
            [](SUBCOIL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(COIL::Name)-1);
            }, "Name of the subcoil")
    ;

    py::class_<COIL>(m, "COIL")
        .def(py::init(&new_Coil), "Create COIL")
        .def_readwrite("Enabled", &COIL::Enabled)
        .def_property("Name", [](COIL& self) {return self.Name;},
            [](COIL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(COIL::Name)-1);
            }, "Name of the coil")
        .def_property_readonly("SubCoils", [](COIL& self) {
            return ObjVecView<SUBCOIL>(self.NumSubCoils, self.SubCoils);},
            py::return_value_policy::reference_internal, "Return vector of SubCoils")
        .def_property("NumSubCoils", [](COIL& self) {return self.NumSubCoils;}, 
            &set_NumSubCoils, "Number of subcoils, setting will erase old subcoils")
        .def_readwrite("CoilCurrent", &COIL::CoilCurrent)
    ;

    py::class_<SUBSHELL>(m, "SUBSHELL")
        .def(py::init(&new_SubShell), "Create SUBSHELL")
        .def_property("Name", [](SUBSHELL& self) {return self.Name;},
            [](SUBSHELL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(SUBSHELL::Name)-1);
            }, "Name of the subshell")
        .def_readwrite("R", &SUBSHELL::X)
        .def_readwrite("Z", &SUBSHELL::Z)
        .def_readwrite("Radius", &SUBSHELL::Radius)
        .def_readwrite("Current", &SUBSHELL::Current)
    ;

    py::class_<SHELL>(m, "SHELL")
        .def(py::init(&new_Shell), "Create SHELL")
        .def_readwrite("Enabled", &SHELL::Enabled)
        .def_property("Name", [](SHELL& self) {return self.Name;},
            [](SHELL& self, std::string name){
                strncpy(self.Name, name.c_str(), sizeof(SHELL::Name)-1);
            }, "Name of the shell")
        .def_readonly("NumSubShells", &SHELL::NumSubShells)
        .def("set_NumSubShells", &set_SHELL_NumSubShells, "Set the number of subshells")
        .def_property_readonly("SubShells", 
            [](SHELL& self) {return ObjVecView<SUBSHELL>(self.NumSubShells, self.SubShells);},
             "Return vector of SubShells")
    ;
}
