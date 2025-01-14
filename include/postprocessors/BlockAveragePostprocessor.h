#pragma once

#include "GeneralPostprocessor.h"
#include "BlockAverageValue.h"
// Forward Declarations
class BlockAveragePostprocessor;

template <>
InputParameters validParams<BlockAveragePostprocessor>();

/**
 *
 */
class BlockAveragePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  BlockAveragePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  unsigned int _block;
  const BlockAverageValue & _block_average_value;
};
