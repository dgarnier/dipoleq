// header classes to make PyDipoleEq a little more readable

#include <pybind11/pybind11.h>

#include "plasma.h"
#include "measurement.h"
#include "tokamak.h"

namespace py = pybind11;


enum class ModelType {
    Std = Plasma_Std,
    IsoNoFlow = Plasma_IsoNoFlow,
    IsoFlow = Plasma_IsoFlow,
    AnisoNoFlow = Plasma_AnisoNoFlow,
    AnisoFlow = Plasma_AnisoFlow,
    DipoleStd = Plasma_DipoleStd,
    DipoleIntStable = Plasma_DipoleIntStable,
    DipoleStablePsiN = Plasma_DipoleStablePsiN
};

enum class MeasType {
    unk = meas_unk,
    bp = meas_bp,
    press = meas_press,
    pperp = meas_pperp,
    ppar = meas_ppar,
    flux = meas_flux,
    saddle = meas_saddle,
    circle = meas_circle,
    coilcur = meas_coilcur,
    plasmacur = meas_plasmacur,
    bt = meas_bt,
    diam = meas_diam,
    bangle = meas_bangle,
    flowt = meas_flowt,
    flowp = meas_flowp,
    ne = meas_ne,
    Te = meas_Te,
    Zeff = meas_Zeff,
    Ti = meas_Ti,
    rho = meas_rho,
    rot = meas_rot,
    ppsix = meas_ppsix,
    bpangle = meas_bpangle,
    pnorm = meas_pnorm,
    J0 = meas_J0
};

enum class CircleType {
    btcos = CircleType_btcos,
    brsin = CircleType_brsin,
    brcos = CircleType_brcos
};

class DMatrixView {
public:
    // matrix is actually nsize+1 x nsize+1
    DMatrixView(std::size_t nsize, double ** data) : m_size(nsize+1), m_data(data), m_dims(2, nsize+1) {};
    double& operator()(std::size_t i, std::size_t j) { return m_data[i][j]; };
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

    static py::object create(std::size_t nsize, double ** data) {
        if (data == NULL)
            return py::none();
        return py::cast(new DMatrixView(nsize, data));
    }
    const std::vector<std::size_t>& Shape() const {
        return m_dims;
    }
    
private:
    std::vector<std::size_t> m_dims;
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
    static py::object create(size_t nsize, int ** data) {
        if (data == NULL)
            return py::none();
        return py::cast(new IMatrixView(nsize, data));
    }
private:
    size_t m_size;
    int ** m_data;
};

class DVectorView {
public:
    // vector size is truthfully nsize
    DVectorView(size_t nsize, double * data) : m_size(nsize), m_data(data) {};
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
    static py::object create(size_t nsize, double * data) {
        if (data == NULL)
            return py::none();
        return py::cast(new DVectorView(nsize, data));
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
        if (i<0 || i>=m_size)
            throw std::out_of_range("Index out of range");
        return m_data[i]; };
    size_t size() {return m_size;};

    void insert(size_t i, T* obj) {
        if (i<0 || i>=m_size)
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

    static py::object create(size_t nsize, T ** data) {
        if (data == NULL)
            return py::none();
        return py::cast(new ObjVecView<T>(nsize, data));
    }
    static py::object create(size_t nsize, T ** data, TOKAMAK * mach, void (*objfree)(T *, TOKAMAK *)) {
        if (data == NULL)
            return py::none();
        return py::cast(new ObjVecView<T>(nsize, data, mach, objfree));
    }

    private:
    TOKAMAK *m_machine = NULL;
    void (*m_free)(T *, TOKAMAK *) = NULL;
    size_t m_size;
    T ** m_data;
};
