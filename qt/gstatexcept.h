#ifndef INCLUDED_GSTATEXCEPT
#define INCLUDED_GSTATEXCEPT



#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif



//namespace pack {



/*!
  \class GStatError
  \brief The GStatError class is a base class for exception classes used
         in the GStat library.

  This base class is useful only for catching *all* errors encountered in the
  GStat lib. Most of the time the derived versions will be used.
*/
//       1         2         3         4         5         6         7         8
class GStatError
{

private:

  //! Message string.
  std::string      d_m;

public:

  //----------------------------------------------------------------------------
  // CREATORS
  //----------------------------------------------------------------------------

  //! Constructor takes a messages string \a m (possibly of type const char *).
                   GStatError          (const std::string &m);

  //! Destructor.
  virtual          ~GStatError         ();

  //----------------------------------------------------------------------------
  // MANIPULATORS
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // ACCESSORS
  //----------------------------------------------------------------------------

  //! Returns the message string.
  const std::string &mesg              () const;

};



//------------------------------------------------------------------------------
// INLINE FUNCIONS
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// FREE OPERATORS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

/*!
  \class GStatNormalError
  \brief The GStatNormalError class is for errors which leave GStat in a
         well defined state.

  If a GStatNormalError occurs GStat can still be used.
*/
//       1         2         3         4         5         6         7         8
class GStatNormalError: public GStatError
{

private:

public:

  //----------------------------------------------------------------------------
  // CREATORS
  //----------------------------------------------------------------------------

  //! Constructor takes a message \a m.
                   GStatNormalError    (const std::string &m);

  //! Destructor.
                   ~GStatNormalError   ();

};

//------------------------------------------------------------------------------

/*!
  \class GStatFatalError
  \brief The GStatFatalError class is for errors which leave GStat in a
         undifined state.

  Use this exception class if a fatal error has occured and GStat can't be
  used anymore.
*/
//       1         2         3         4         5         6         7         8
class GStatFatalError: public GStatError
{

private:

public:

  //----------------------------------------------------------------------------
  // CREATORS
  //----------------------------------------------------------------------------

  //! Constructor takes a message \a m.
                   GStatFatalError     (const std::string &m);

  //! Destructor.
                   ~GStatFatalError    ();

};

//} // namespace pack

#endif
