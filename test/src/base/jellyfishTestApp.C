//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "jellyfishTestApp.h"
#include "jellyfishApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
jellyfishTestApp::validParams()
{
  InputParameters params = jellyfishApp::validParams();
  return params;
}

jellyfishTestApp::jellyfishTestApp(InputParameters parameters) : MooseApp(parameters)
{
  jellyfishTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

jellyfishTestApp::~jellyfishTestApp() {}

void
jellyfishTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  jellyfishApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"jellyfishTestApp"});
    Registry::registerActionsTo(af, {"jellyfishTestApp"});
  }
}

void
jellyfishTestApp::registerApps()
{
  registerApp(jellyfishApp);
  registerApp(jellyfishTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
jellyfishTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  jellyfishTestApp::registerAll(f, af, s);
}
extern "C" void
jellyfishTestApp__registerApps()
{
  jellyfishTestApp::registerApps();
}
