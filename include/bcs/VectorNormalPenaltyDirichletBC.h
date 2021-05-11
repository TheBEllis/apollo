//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorIntegratedBC.h"

class VectorNormalPenaltyDirichletBC : public VectorIntegratedBC
{
public:
  static InputParameters validParams();

  VectorNormalPenaltyDirichletBC(const InputParameters & parameters);

protected:
  Real _penalty;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  // Holds the solution curl at the current quadrature points
  const Function * const _function;
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;

};