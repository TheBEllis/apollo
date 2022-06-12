#pragma once

#include "mfem.hpp"

class MFEMMesh : public mfem::Mesh {
 public:
  MFEMMesh(int num_elem, std::vector<double> coordx, std::vector<double> coordy, std::vector<double> coordz, 
          std::map<int, int> cubitToMFEMVertMap, std::vector<int> uniqueVertexID,
          int libmesh_element_type, int libmesh_face_type, int** elem_blk, int num_el_blk,
          unsigned int num_node_per_el, size_t* num_el_in_blk, int num_element_linear_nodes, int num_face_nodes,
          int num_face_linear_nodes, int num_side_sets, std::vector<int> num_sides_in_ss, int** ss_node_id, int* ebprop,
          int * ssprop /*,3*/);

  MFEMMesh(std::string afilename, int generate_edges = 0, int refine = 1,
           bool fix_orientation = true);

  enum CubitFaceType {
    FACE_EDGE2,
    FACE_EDGE3,
    FACE_TRI3,
    FACE_TRI6,
    FACE_QUAD4,
    FACE_QUAD9
  };

  enum CubitElementType {
    ELEMENT_TRI3,
    ELEMENT_TRI6,
    ELEMENT_QUAD4,
    ELEMENT_QUAD9,
    ELEMENT_TET4,
    ELEMENT_TET10,
    ELEMENT_HEX8,
    ELEMENT_HEX27
  };

};
