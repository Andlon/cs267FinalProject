#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GSTATEXCEPT
#include "gstatexcept.h"
#endif



GStatError::GStatError(const std::string &m)

  : d_m(m)

{
}



GStatError::~GStatError()
{
}



const std::string &GStatError::mesg() const
{
  return d_m;
}

//------------------------------------------------------------------------------

GStatNormalError::GStatNormalError(const std::string &m)

  : GStatError(m)

{
}



GStatNormalError::~GStatNormalError()
{
}

//------------------------------------------------------------------------------

GStatFatalError::GStatFatalError(const std::string &m)

  : GStatError(m)

{
}



GStatFatalError::~GStatFatalError()
{
}

