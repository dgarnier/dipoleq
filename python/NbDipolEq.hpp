// NbDipolEq.hpp - nanobind helper classes for DipolEq
// New nanobind binding; see PyDipolEq.hpp for the pybind11 equivalent (kept for reference).
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
// #include <nanobind/stl/optional.h>
// #include <optional>
#include <stdexcept>

#include "plasma.h"
#include "measurement.h"
#include "tokamak.h"

namespace nb = nanobind;

// ── Enum types (same underlying values as PyDipolEq.hpp) ──────────────────────

enum class ModelType {
    Std              = Plasma_Std,
    IsoNoFlow        = Plasma_IsoNoFlow,
    IsoFlow          = Plasma_IsoFlow,
    AnisoNoFlow      = Plasma_AnisoNoFlow,
    AnisoFlow        = Plasma_AnisoFlow,
    DipoleStd        = Plasma_DipoleStd,
    DipoleIntStable  = Plasma_DipoleIntStable,
    DipoleStablePsiN = Plasma_DipoleStablePsiN
};

enum class MeasType {
    unk       = meas_unk,
    bp        = meas_bp,
    press     = meas_press,
    pperp     = meas_pperp,
    ppar      = meas_ppar,
    flux      = meas_flux,
    saddle    = meas_saddle,
    circle    = meas_circle,
    coilcur   = meas_coilcur,
    plasmacur = meas_plasmacur,
    bt        = meas_bt,
    diam      = meas_diam,
    bangle    = meas_bangle,
    flowt     = meas_flowt,
    flowp     = meas_flowp,
    ne        = meas_ne,
    Te        = meas_Te,
    Zeff      = meas_Zeff,
    Ti        = meas_Ti,
    rho       = meas_rho,
    rot       = meas_rot,
    ppsix     = meas_ppsix,
    bpangle   = meas_bpangle,
    pnorm     = meas_pnorm,
    J0        = meas_J0
};

enum class CircleType {
    btcos = CircleType_btcos,
    brsin = CircleType_brsin,
    brcos = CircleType_brcos
};

// ── Array view helpers ────────────────────────────────────────────────────────
//
// These return numpy-compatible ndarray views into C-managed memory.
// The data layout matches the pybind11 version (nrutil "Fortran-order" strides:
// first index varies fastest).
//
// IMPORTANT: The returned array is only valid while the parent C object is alive.
// The caller is responsible for lifetime management.

using dmatrix_v = nb::ndarray<double, nb::numpy, nb::ndim<2>, nb::f_contig>;
using dvector_v = nb::ndarray<double, nb::numpy, nb::ndim<1>, nb::f_contig>;
using imatrix_v = nb::ndarray<int, nb::numpy, nb::ndim<2>, nb::f_contig>;


inline nb::object make_dmatrix(size_t nsize, double **data) {
    if (!data) return nb::none();
    size_t n = nsize + 1;
    return nb::cast(dmatrix_v(data[0], {n, n}),
        nb::rv_policy::reference_internal);
}

inline nb::object make_dvector(size_t n, double *data) {
    if (!data) return nb::none();
    return nb::cast(dvector_v(data, {n}),
        nb::rv_policy::reference_internal);
}

inline nb::object make_imatrix(size_t nsize, int **data) {
    if (!data) return nb::none();
    size_t n = nsize + 1;
    return nb::cast(imatrix_v(data[0], {n, n}),
        nb::rv_policy::reference_internal);
}

// ── ObjVecView ────────────────────────────────────────────────────────────────
//
// A non-owning view over an array of C object pointers.
// Provides sequence-like access; the underlying pointers are owned by TOKAMAK.

template <typename T>
class ObjVecView {
public:
    ObjVecView() = default;

    ObjVecView(size_t n, T **data)
        : m_size(n), m_data(data) {}

    ObjVecView(size_t n, T **data, TOKAMAK *mach, void (*objfree)(T *, TOKAMAK *))
        : m_size(n), m_data(data), m_machine(mach), m_free(objfree) {}

    T *&operator[](size_t i) {
        if (i >= m_size) throw std::out_of_range("Index out of range");
        return m_data[i];
    }

    size_t size() const { return m_size; }

    // Replace element i (freeing the old one if a deleter was provided).
    void insert(size_t i, T *obj) {
        if (i >= m_size) throw std::out_of_range("Index out of range");
        if (m_data[i]) {
            if (m_free) m_free(m_data[i], m_machine);
            else        free(m_data[i]);
        }
        m_data[i] = obj;
    }
    // give the start and end of array for the iterator
    auto begin() const { return &m_data[0]; }
    auto end() const { return &m_data[m_size]; }

    TOKAMAK *m_machine = nullptr;
private:
    void (*m_free)(T *, TOKAMAK *) = nullptr;
    size_t  m_size = 0;
    T     **m_data = nullptr;
};
