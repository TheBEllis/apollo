#pragma once
#include "MFEMMesh.h"
#include "libmesh/face_quad4.h"

// Constructor to create an MFEM mesh from VTK data structures. These data
// structures are obtained by the methods found in MFEMproblem
MFEMMesh::MFEMMesh(int num_elem, std::vector<double> coordx, std::vector<double> coordy, std::vector<double> coordz, 
      std::map<int, int> cubitToMFEMVertMap, std::vector<int> uniqueVertexID,
      int libmesh_element_type, int libmesh_face_type, int** elem_blk, int num_el_blk,
      unsigned int num_node_per_el, size_t* num_el_in_blk, int num_element_linear_nodes, int num_face_nodes,
      int num_face_linear_nodes, int num_side_sets, std::vector<int> num_sides_in_ss, int** ss_node_id, int* ebprop,
      int * ssprop /*,3*/) {
  SetEmpty();

  NumOfElements = num_elem;
  elements.SetSize(num_elem);
  int elcount = 0;
  int renumberedVertID[8];
  for (int iblk = 0; iblk < (int) num_el_blk; iblk++)
  {
    int NumNodePerEl = num_node_per_el;
    for (int i = 0; i < (int) num_el_in_blk[iblk]; i++)
    {
        for (int j = 0; j < num_element_linear_nodes; j++)
        {
          renumberedVertID[j] =
              cubitToMFEMVertMap[elem_blk[iblk][i*NumNodePerEl+j]]-1;
        }

        switch (cubit_element_type)
        {
          case (ELEMENT_TRI3):
          case (ELEMENT_TRI6):
          {
              elements[elcount] = new Triangle(renumberedVertID,ebprop[iblk]);
              break;
          }
          case (ELEMENT_QUAD4):
          case (ELEMENT_QUAD9):
          {
              elements[elcount] = new Quadrilateral(renumberedVertID,ebprop[iblk]);
              break;
          }
          case (ELEMENT_TET4):
          case (ELEMENT_TET10):
          {
#ifdef MFEM_USE_MEMALLOC
              elements[elcount] = TetMemory.Alloc();
              elements[elcount]->SetVertices(renumberedVertID);
              elements[elcount]->SetAttribute(ebprop[iblk]);
#else
              elements[elcount] = new Tetrahedron(renumberedVertID,
                                                  ebprop[iblk]);
#endif
              break;
          }
          case (ELEMENT_HEX8):
          case (ELEMENT_HEX27):
          {
              elements[elcount] = new Hexahedron(renumberedVertID,ebprop[iblk]);
              break;
          }
        }
        elcount++;
    }
  }

  // load up the boundary elements

  NumOfBdrElements = 0;
  for (int iss = 0; iss < (int) num_side_sets; iss++)
  {
    NumOfBdrElements += num_side_in_ss[iss];
  }
  boundary.SetSize(NumOfBdrElements);
  int sidecount = 0;
  for (int iss = 0; iss < (int) num_side_sets; iss++)
  {
    for (int i = 0; i < (int) num_side_in_ss[iss]; i++)
    {
        for (int j = 0; j < num_face_linear_nodes; j++)
        {
          renumberedVertID[j] =
              cubitToMFEMVertMap[ss_node_id[iss][i*num_face_nodes+j]] - 1;
        }
        switch (cubit_face_type)
        {
          case (FACE_EDGE2):
          case (FACE_EDGE3):
          {
              boundary[sidecount] = new Segment(renumberedVertID,ssprop[iss]);
              break;
          }
          case (FACE_TRI3):
          case (FACE_TRI6):
          {
              boundary[sidecount] = new Triangle(renumberedVertID,ssprop[iss]);
              break;
          }
          case (FACE_QUAD4):
          case (FACE_QUAD9):
          {
              boundary[sidecount] = new Quadrilateral(renumberedVertID,ssprop[iss]);
              break;
          }
        }
        sidecount++;
    }
  }
}

// Constructor to create an MFEM mesh from a file, currently just used for
// testing purposes
MFEMMesh::MFEMMesh(std::string cppfilename, int generate_edges, int refine,
                   bool fix_orientation) {
  std::cout << cppfilename << std::endl;
  const char *filename = cppfilename.c_str();
  SetEmpty();

  mfem::named_ifgzstream imesh(filename);
  if (!imesh) {
    MFEM_ABORT("Mesh file not found: " << filename << '\n');
  } else {
    Load(imesh, generate_edges, refine, fix_orientation);
  }
}

