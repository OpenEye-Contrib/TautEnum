//
// file MolDiffViewer2.H
// David Cosgrove
// AstraZeneca
// 24th December 2013

// This is the main class for program mol_diff_viewer2.

#ifndef MOLDIFFVIEWER2_H
#define MOLDIFFVIEWER2_H
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
class QGridLayout;
class QLabel;
class QLineEdit;
class QSlider;
class QTextEdit;

namespace DACLIB {
class QTMolDisplay2D;
}

// ****************************************************************************

class MolDiffViewer2 : public QMainWindow {

  Q_OBJECT

public :

  MolDiffViewer2( int argc , char **argv );

private :

  QAction *file_quit_ , *next_left_smaller_ , *write_diffs_;
  QSlider *mol_slider_;
  QLabel *file1_label_ , *file2_label_;
  std::vector<DACLIB::QTMolDisplay2D *> mol_disp_ls_ , mol_disp_rs_;
  QGridLayout *left_grid_ , *right_grid_;
  QTextEdit *smiles_1_ , *smiles_2_;
  std::vector<std::pair<std::string,std::vector<std::string> > > file1_mols_ , file2_mols_;

  std::vector<int> diffs_;

  void build_widget();
  void build_actions();
  void build_menubar();

  void read_smiles_file( const std::string &filename ,
                         std::vector<std::pair<std::string,std::vector<std::string> > > &mols );
  void find_differences();

  void show_smiles( std::pair<std::string,std::vector<std::string> > &mols ,
                    std::vector<DACLIB::QTMolDisplay2D *> &mol_disps ,
                    QGridLayout *grid , QTextEdit *other_smiles );

private slots :

  void slot_quit();
  void slot_slider_changed();
  void slot_next_left_smaller();
  void slot_write_diffs();

};

#endif // MOLDIFFVIEWER_H

#endif // MOLDIFFVIEWER2_H
