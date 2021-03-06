#include "jellyfishApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
jellyfishApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

jellyfishApp::jellyfishApp(InputParameters parameters) : MooseApp(parameters)
{
  jellyfishApp::registerAll(_factory, _action_factory, _syntax);
}

jellyfishApp::~jellyfishApp() {}

void
jellyfishApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"jellyfishApp"});
  Registry::registerActionsTo(af, {"jellyfishApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
jellyfishApp::registerApps()
{
  registerApp(jellyfishApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
jellyfishApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  jellyfishApp::registerAll(f, af, s);
}
extern "C" void
jellyfishApp__registerApps()
{
  jellyfishApp::registerApps();
}
