#include "pbxApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
pbxApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

pbxApp::pbxApp(InputParameters parameters) : MooseApp(parameters)
{
  pbxApp::registerAll(_factory, _action_factory, _syntax);
}

pbxApp::~pbxApp() {}

void
pbxApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"pbxApp"});
  Registry::registerActionsTo(af, {"pbxApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
pbxApp::registerApps()
{
  registerApp(pbxApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
pbxApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  pbxApp::registerAll(f, af, s);
}
extern "C" void
pbxApp__registerApps()
{
  pbxApp::registerApps();
}
