//
// file TautEnumCallableThreaded.H
// David Cosgrove
// AstraZeneca
// 23rd February 2012
//
// This is a concrete class of TautEnumCallable which does tautomer enumeration
// via threads. It's necessary because the OEChem threaded I/O is not
// in the same class as the non-threaded, so you can't just use a pointer to
// the stream objects and pass the appropriate one as required. Seems like a
// design flaw to me, but what do I really know about threading?

#ifndef TAUTENUMCALLABLETHREADED_H
#define TAUTENUMCALLABLETHREADED_H

#include "TautEnumCallableBase.H"

#include <oechem.h>

// ****************************************************************************

class TautEnumCallableThreaded : public TautEnumCallableBase {

public :

  TautEnumCallableThreaded( OEChem::oemolithread *in , OEChem::oemolothread *out ,
                            const TautEnumSettings &settings ) :
    TautEnumCallableBase( settings ) ,
    in_( in ) , out_( out ) {}

  TautEnumCallableThreaded( const TautEnumCallableThreaded &rhs ) :
    TautEnumCallableBase( rhs.tes_ ) ,
    in_( rhs.in_ ) , out_( rhs.out_ ) {}

  ~TautEnumCallableThreaded() {}

private :

  // oemol[io]thread are mutexed properly so there should be no problems
  // using multiple threads to read from and write to the same stream
  OEChem::oemolithread *in_;
  OEChem::oemolothread *out_;

  bool read_next_molecule( OEChem::OEMolBase &mol ) {
    return OEChem::OEReadMolecule( *in_ , mol );
  }
  void write_molecule( OEChem::OEMolBase &mol ) {
    OEChem::OEWriteMolecule( *out_ , mol );
  }

};

#endif // TAUTENUMCALLABLETHREADED_H
