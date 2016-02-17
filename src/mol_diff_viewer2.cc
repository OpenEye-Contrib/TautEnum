//
// file mol_diff_viewer2.cc
// David Cosgrove
// AstraZeneca
// 22nd November 2013
//
// A quick and dirty program to show differences in molecules from two files.
// It takes 2 files, with molecules assumed to be in the same order by size,
// and shows all molecules of the same name where there are differences between
// the SMILES strings. Multiple molecules of the same name are expected.

#include "MolDiffViewer2.H"

#include <QApplication>

// ****************************************************************************
int main( int argc , char **argv ) {

  QApplication a( argc , argv );

  MolDiffViewer2 mdf2( argc , argv );
  mdf2.setGeometry( 100 , 100 , 900 , 500 );
  mdf2.show();

  return a.exec();

}

