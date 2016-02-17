//
// file mol_diff_viewer.cc
// David Cosgrove
// AstraZeneca
// 22nd November 2013
//
// A quick and dirty program to show differences in molecules from two files.
// It will show side-by-side all molecule pairs, one from each of two files
// that have the name but different SMILES strings/structures.  The purpose
// being to allow easy comparison of the differences between output from
// two different structure standardisation programs, by focussing on the
// cases where the results are different.  To that end, it assumes that
// the two input files have the same molecules in the same order.
// Takes 2 or 4 command-line args.  The first 2 are the two files as
// described.  If 4 are given, the 2nd pair correspond to the 1st 2,
// but are an enumerated set of tautomers for the same molecules,
// with all SMILES strings for the different molecules being shown,
// so one can decide which standardisation program did the better
// enumeration job.

#include "MolDiffViewer.H"

#include <QApplication>

// ****************************************************************************
int main( int argc , char **argv ) {

  QApplication a( argc , argv );

  MolDiffViewer mdf( argc , argv );
  mdf.setGeometry( 100 , 100 , 900 , 500 );
  mdf.show();

  return a.exec();

}
