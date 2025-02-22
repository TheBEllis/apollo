//* Base class for solution of (real) Maxwell equations in time domain.

#include "MaxwellBase.h"
#include "Function.h"
#include "Assembly.h"

registerMooseObject("ApolloApp", MaxwellBase);

InputParameters
MaxwellBase::validParams()
{
  InputParameters params = VectorTimeKernel::validParams();
  params.addClassDescription("This class computes various components of the"
                             "Maxwell equations which can then be assembled "
                             "together in child classes.");
  params.addParam<bool>("quasistationary", true, "Whether the term proportional to the second order time derivative of the field should be neglected."
  "If false, a time integrator such as NewmarkBeta that calculates second order time derivatives must be used.");
  params.addParam<bool>("gauge_penalty",
                        false,
                        "Whether a gauge penalty term should be included for regularisation.");

  return params;
}

MaxwellBase::MaxwellBase(const InputParameters & parameters)
  : VectorTimeKernel(parameters),
    _curl_phi(_assembly.curlPhi(_var)),
    _curl_test(_var.curlPhi()),
    _curl_u(_is_implicit ? _var.curlSln() : _var.curlSlnOld()),
    _u_dotdot(getParam<bool>("quasistationary") ?   _vector_zero :_var.uDotDot()),
    _du_dotdot_du(getParam<bool>("quasistationary") ?  _zero : _var.duDotDotDu()),
    _quasistationary(getParam<bool>("quasistationary")),
    _gauge_penalty(getParam<bool>("gauge_penalty"))
{
}

Real
MaxwellBase::curlCurlTerm()
{
  return _curl_u[_qp] * _curl_test[_i][_qp];
}

Real
MaxwellBase::dCurlCurlDU()
{
  return _curl_phi[_j][_qp] * _curl_test[_i][_qp];
}

Real
MaxwellBase::gaugePenaltyTerm()
{
  return _grad_u[_qp].tr() * _grad_test[_i][_qp].tr();
}

Real
MaxwellBase::dGaugePenaltyDU()
{
  return _grad_phi[_j][_qp].tr() * _grad_test[_i][_qp].tr();
}

Real
MaxwellBase::steadyStateTerm()
{
  return curlCurlTerm() + (_gauge_penalty ? gaugePenaltyTerm() : 0);
}

Real
MaxwellBase::dSteadyStateTermDU()
{
  return dCurlCurlDU() + (_gauge_penalty ? dGaugePenaltyDU() : 0);
}

Real
MaxwellBase::firstOrderTimeDerivTerm()
{
  return _test[_i][_qp] * _u_dot[_qp];
}

Real
MaxwellBase::dFirstOrderTimeDerivDU()
{
  return _test[_i][_qp] * _du_dot_du[_qp] * _phi[_j][_qp];
}

Real
MaxwellBase::secondOrderTimeDerivTerm()
{
  return _test[_i][_qp] * _u_dotdot[_qp];
}

Real
MaxwellBase::dSecondOrderTimeDerivDU()
{
  return _test[_i][_qp] *  _du_dotdot_du[_qp] * _phi[_j][_qp];
}

Real
MaxwellBase::computeQpResidual() {return 0;}

Real
MaxwellBase::computeQpJacobian(){return 0;}