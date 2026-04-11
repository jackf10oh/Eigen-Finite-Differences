// PyFds.cpp
//
// pybind11 bindsing for LinOps mesh + discretization
//
// JAF 1/25/2026 

#include<pybind11/pybind11.h> 
#include<pybind11/stl.h>
#include<pybind11/functional.h> 
#include<pybind11/eigen.h>

#include<LinOps/All.hpp>
#include "Wave2D.hpp" 

namespace py = pybind11; 

PYBIND11_MODULE(PyFds, m)
{
  // Mesh1D ================================================================================ 
  py::class_<LinOps::Mesh1D, std::shared_ptr<LinOps::Mesh1D>>(m, "Mesh1D")
    .def(py::init<double,double,std::size_t>(), 
           py::arg("start")=0.0, 
           py::arg("stop")=1.0, 
           py::arg("nsteps")=11
    )
    .def(py::init<std::vector<double>>(), 
           py::arg("vals")
    )
    .def("at",[](const LinOps::Mesh1D& self, std::size_t i){return self.at(i);}, py::arg("idx"), "values of mesh at an index") 
    .def("size", &LinOps::Mesh1D::size, "number of entries in a Mesh1D")   
    .def("__len__", &LinOps::Mesh1D::size) // combines with __getitem__ to make mesh1d iterable in python... 
    .def("__getitem__", [](const LinOps::Mesh1D& mesh, size_t i) {
        if (i >= mesh.size()) throw py::index_error();
        return mesh[i];  // return the coordinate at index i
    });

  // MeshXD ================================================================================ 
  py::class_<LinOps::MeshXD, std::shared_ptr<LinOps::MeshXD>>(m,"MeshXD")
    .def(py::init<double,double,std::size_t,std::size_t>(), 
           py::arg("start")=0.0, 
           py::arg("stop")=1.0, 
           py::arg("nsteps")=11, 
           py::arg("ndims")=2
    )
    .def(py::init<std::vector<std::pair<double,double>>, std::vector<std::size_t>>(), 
           py::arg("endpoints_list"), 
           py::arg("nsteps_list")
    )
    .def(py::init<std::shared_ptr<const LinOps::Mesh1D>>(),
      py::arg("mesh1d")
    )
    .def(py::init<std::vector<std::shared_ptr<const LinOps::Mesh1D>>>(),
           py::arg("mesh1d_list")
    )
    .def("getMesh1D", 
           [](LinOps::MeshXD& self, std::size_t i){ return self.getMesh1DSafe(i); },
           py::arg("dim")=1)
    .def("numDims", 
           &LinOps::MeshXD::numDims)
    .def("sizeOfDim", 
           &LinOps::MeshXD::sizeOfDim,
           py::arg("dim")=1)
    .def("sizesProduct", 
           &LinOps::MeshXD::sizesProduct)
    .def("sizesMiddleProduct", 
           &LinOps::MeshXD::sizesMiddleProduct,
           py::arg("start"),
           py::arg("stop"));
  
  // LinOps::Vector1D ================================================================================ 
  py::class_<LinOps::Vector1D, std::unique_ptr<LinOps::Vector1D>>(m,"Vector1D")
    .def(py::init<std::size_t>(), 
           py::arg("size")=0
    )
    .def(py::init<std::shared_ptr<const LinOps::Mesh1D>>(), 
           py::arg("mesh")
    )
    .def(py::init<Eigen::VectorXd&&>(),
           py::arg("arr")
    )
    .def(py::init<const Eigen::VectorXd&>(),
           py::arg("arr")
    )
    .def("values",
            [](const LinOps::Vector1D& self){return self.values();}, 
            py::return_value_policy::reference_internal
    )
    .def("size", 
            &LinOps::Vector1D::size, 
            "size of current Vector1D"
    )
    .def("resize",
            &LinOps::Vector1D::resize, 
            py::arg("mesh"),"resize Vector1D to fit on a Mesh1D"
    ); 

  // LinOps::VectorXD ================================================================================ 
  py::class_<LinOps::VectorXD, std::unique_ptr<LinOps::VectorXD>>(m,"VectorXD")
    .def(py::init<std::size_t>(), 
           py::arg("size")=0
    )
    .def(py::init<std::shared_ptr<const LinOps::MeshXD>>(), 
           py::arg("mesh")
    )
    .def(py::init<Eigen::VectorXd&&>(),
           py::arg("arr")
    )
    .def(py::init<const Eigen::VectorXd&>(),
           py::arg("arr")
    )
    .def("values",
            [](const LinOps::VectorXD& self){return self.values();}, 
            py::return_value_policy::reference_internal
    )
    .def("sizesProduct", 
            &LinOps::VectorXD::sizesProduct, 
            "product of each dimensions size"
    )
    .def("numDims", 
            &LinOps::VectorXD::numDims, 
            "number of dimensions in VectorXD"
    )
    .def("sizeOfDim", 
            &LinOps::VectorXD::sizeOfDim,
            py::arg("dim")=0, 
            "size of specific dimension in VectorXD"
    )
    .def("sizesMiddleProduct", 
            &LinOps::VectorXD::sizesMiddleProduct,
            py::arg("start"), 
            py::arg("stop"), 
            "product of sizes of dimensions in [start,stop)" 
    )
    .def("resize",
            &LinOps::VectorXD::resize, 
            py::arg("mesh"),"resize VectorXD to fit on a MeshXD"
    ); 

  // Wave2D ================================================================================ 
  py::class_<Wave2D>(m, "Wave2D")
    .def(py::init<>())
    .def("SetDomain", [](Wave2D& self, LinOps::SharedConstMeshXD m) -> Wave2D& 
      {auto a = self.Args(); a.domain_mesh_ptr=m; self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new domain mesh inside of solver") 
    .def("SetTime", [](Wave2D& self, LinOps::SharedConstMesh1D m) -> Wave2D& 
      {auto a = self.Args(); a.time_mesh_ptr=m; self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new time mesh inside of solver") 
    .def("SetIC", [](Wave2D& self, std::vector<Eigen::VectorXd> v) -> Wave2D& 
      {auto a = self.Args(); a.ICs=std::move(v); self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new initial condition inside of solver") 
    .def("SetHeight",[](Wave2D& self, double h) -> Wave2D& 
      { self.set_bump_height(h); self.Reset(); return self;}, py::arg("h")=1.0, "Sets new height for oscilation force term at origin")
    .def("SetDamping", [](Wave2D& self, double d) -> Wave2D& 
      {self.set_damping(d); self.Reset(); return self;}, py::arg("d")=2.0, "Sets a new damping rate at boundaries of the domain")
    .def("StoredData", &Wave2D::StoredData, "returns a vector of solutions. solutions are flattened in dimensional order")
    .def("Compute", &Wave2D::FillVals, py::arg("max_iters")=20, "Computes solution at each entry in time. available in StoredData")
    .def("SolAt", [](Wave2D& self, double t, double x, double y){ return self.SolAt(t,x,y); }, py::arg("t"), py::arg("x"), py::arg("y"), "Returns value of solution at time t at coords (x,y)");    
}

