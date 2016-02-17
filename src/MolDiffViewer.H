//
// file MolDiffViewer.H
// David Cosgrove
// AstraZeneca
// 22nd November 2013

// This is the main class for program mol_diff_viewer.

#ifndef MOLDIFFVIEWER_H
#define MOLDIFFVIEWER_H

// ****************************************************************************

#include <string>
#include <vector>

#include <QMainWindow>

class QAction;
class QLabel;
class QLineEdit;
class QSlider;
class QTextEdit;

namespace DACLIB {
class QTMolDisplay2D;
}

// ****************************************************************************

class MolDiffViewer : public QMainWindow {

  Q_OBJECT

public :

  MolDiffViewer( int argc , char **argv );

private :

  QAction *file_quit_;
  QSlider *mol_slider_;
  QLabel *file1_label_ , *file2_label_;
  DACLIB::QTMolDisplay2D *mol_disp_l_ , *mol_disp_r_;
  QLineEdit *smiles_1_ , *smiles_2_;
  QTextEdit *other_smiles_1_ , *other_smiles_2_;
  std::vector<std::pair<std::string,std::string> > file1_mols_ , file2_mols_;
  std::vector<std::pair<std::string,std::string> > file3_mols_ , file4_mols_;

  std::vector<int> diffs_;

  void build_widget();
  void build_actions();
  void build_menubar();

  void read_smiles_file( const std::string &filename ,
                         std::vector<std::pair<std::string,std::string> > &mols );
  void find_differences();
  void show_other_smiles( const std::string &mol_name ,
                          const std::vector<std::pair<std::string,std::string> > &mols ,
                          QTextEdit *other_smiles );

private slots :

  void slot_quit();
  void slot_slider_changed();

};

#endif // MOLDIFFVIEWER_H
