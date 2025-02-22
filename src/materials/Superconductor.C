#include "Superconductor.h"
#include "Function.h"

registerMooseObject("ApolloApp", Superconductor);

InputParameters
Superconductor::validParams()
{
  InputParameters params = Material::validParams();

  params.addParam<Real>(
      "critical_electric_field",
      1.0,
      "The critical electric field ($E_c$) of the superconductor");
  params.addParam<Real>(
      "critical_current_density",
      1.0,
      "The critical current density ($J_c$) of the superconductor");
  params.addParam<Real>(
      "permeability", 1.0, "The permeability ($\\mu$) of the conductor. Defaults to 1");
  params.addParam<Real>("nonlinearity_parameter",
                        1.0,
                        "The nonlinearity parameter ($n$) of the superconductor. Defaults to 1, "
                        "representing a linear resistivity term");
  params.addCoupledVar("magnetic_field", "The magnetic field ($H$) as a function of position.");

  return params;
}

Superconductor::Superconductor(const InputParameters & parameters)
  : Material(parameters),
    // DerivativeMaterialInterface<Material>(parameters),
    // Get the parameters from the input file
    _input_n(getParam<Real>("nonlinearity_parameter")),
    _input_jc(getParam<Real>("critical_current_density")),
    _input_ec(getParam<Real>("critical_electric_field")),
    _input_permeability(getParam<Real>("permeability")),
    // Declare material properties
    _n(declareProperty<Real>("nonlinearity_parameter")),
    _jc(declareProperty<Real>("critical_current_density")),
    _ec(declareProperty<Real>("critical_electric_field")),
    // _H(coupledValue("magnetic_field")),
    // _H_var(coupled("magnetic_field")),
    // _H_name(getVar("magnetic_field", 0)->name()),
    // Declare two material properties by getting a reference from the MOOSE Material system
    _permeability(declareProperty<Real>("permeability")),
    _resistivity(declareProperty<Real>("resistivity")),
    _drhodj(declareProperty<Real>("drhodj")),
    // _drdH(declarePropertyDerivative<Real>("resistivity", _H_name)),
    // _d2rdH2(declarePropertyDerivative<Real>("resistivity", _H_name, _H_name))

    // Declare two material properties by getting a reference from the MOOSE Material system
    // _h(isCoupled("magnetic_field"),
    _j(coupledCurl("magnetic_field"))
{
}

void
Superconductor::computeQpProperties()
{
  _ec[_qp] = _input_ec;
  _jc[_qp] = _input_jc;
  _n[_qp] = _input_n;
  _permeability[_qp] = _input_permeability;
  _resistivity[_qp] = (_ec[_qp] / _jc[_qp]) * pow((_j[_qp].norm() / _jc[_qp]), _n[_qp] - 1);
  _drhodj[_qp] = (_ec[_qp] / _jc[_qp]) * _n[_qp] * (pow((_j[_qp].norm() / _jc[_qp]), _n[_qp] - 1));
}
