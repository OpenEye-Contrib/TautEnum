//
// file MolDiffViewer.cc
// David Cosgrove
// AstraZeneca
// 22nd November 2013

#include "MolDiffViewer.H"

#include "QTMolDisplay2D.H"

#include <QAction>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QSlider>
#include <QStatusBar>
#include <QString>
#include <QTextEdit>

#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>

#include <oechem.h>

using namespace OEChem;
using namespace std;

// ****************************************************************************
MolDiffViewer::MolDiffViewer( int argc , char **argv ) : QMainWindow() {

  build_widget();

  read_smiles_file( argv[1] , file1_mols_ );
  file1_label_->setText( argv[1] );
  read_smiles_file( argv[2] , file2_mols_ );
  file2_label_->setText( argv[2] );

  if( argc > 3 ) {
    read_smiles_file( argv[3] , file3_mols_ );
    read_smiles_file( argv[4] , file4_mols_ );
  } else {
    other_smiles_1_->hide();
    other_smiles_2_->hide();
  }

  find_differences();
  if( !diffs_.empty() ) {
    mol_slider_->setMaximum( diffs_.size() - 1 );
    mol_slider_->setMinimum( 0 );
  } else {
    mol_slider_->setEnabled( false );
  }

  slot_slider_changed();

}

// ****************************************************************************
void MolDiffViewer::build_widget() {

  build_actions();
  build_menubar();

  mol_disp_l_ = new DACLIB::QTMolDisplay2D;
  mol_disp_r_ = new DACLIB::QTMolDisplay2D;
  mol_slider_ = new QSlider;
  mol_slider_->setPageStep( 1 );
  connect( mol_slider_ , SIGNAL( valueChanged(int) ) ,
           this , SLOT( slot_slider_changed() ) );

  QVBoxLayout *vbox1 = new QVBoxLayout;
  vbox1->addWidget( mol_disp_l_ , 1 );
  smiles_1_ = new QLineEdit;
  vbox1->addWidget( smiles_1_ );
  other_smiles_1_ = new QTextEdit;
  vbox1->addWidget( other_smiles_1_ );
  file1_label_ = new QLabel;
  vbox1->addWidget( file1_label_ );

  QVBoxLayout *vbox2 = new QVBoxLayout;
  vbox2->addWidget( mol_disp_r_ , 1 );
  smiles_2_ = new QLineEdit;
  vbox2->addWidget( smiles_2_ );
  other_smiles_2_ = new QTextEdit;
  vbox2->addWidget( other_smiles_2_ );
  file2_label_ = new QLabel;
  vbox2->addWidget( file2_label_ );

  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addLayout( vbox1 );
  hbox->addWidget( mol_slider_ );
  hbox->addLayout( vbox2 );

  QWidget *central_wid = new QWidget;
  central_wid->setLayout( hbox );

  setCentralWidget( central_wid );

}

// ****************************************************************************
void MolDiffViewer::build_actions() {

  file_quit_ = new QAction( "Quit" , this );
  file_quit_->setShortcut( QString( "Ctrl+Q" ) );
  connect( file_quit_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_quit() ) );

}

// ****************************************************************************
void MolDiffViewer::build_menubar() {

  QMenu *file_menu = menuBar()->addMenu( "File" );
  file_menu->addAction( file_quit_ );

}

// ****************************************************************************
void MolDiffViewer::read_smiles_file( const string &filename ,
                                      vector<pair<string,string> > &mols ) {

  using namespace boost;

  cout << "Reading " << filename << endl;

  ifstream file( filename.c_str() , ios_base::in | ios_base::binary );
  iostreams::filtering_streambuf<boost::iostreams::input> ins;
  if( filename.length() > 7 && filename.substr( filename.length() - 7 ) == string( ".smi.gz" ) ) {
    ins.push( iostreams::gzip_decompressor() );
  }
  ins.push( file );

  istreambuf_iterator<char> in( &ins ) , eos;

  vector<char> nextline;

  int mol_num = 1;
  while( in != eos ) {
    nextline.clear();
    while( in != eos && *in != '\n' ) {
      nextline.push_back( *in );
      ++in;
    }
    vector<string> split_line;
    split( split_line , nextline , is_any_of( " \t," ) , boost::algorithm::token_compress_on );
    // take the molecule name as anything after the first field, so that spaces can be in
    // the name
    string mol_name;
    if( split_line.size() > 1 ) {
      mol_name = string( nextline.begin() + split_line[0].length() + 1 ,
                         nextline.end() );
      boost::trim( mol_name );
    } else {
      mol_name = string( "Mol_" ) + lexical_cast<string>( mol_num );
    }
    boost::trim( split_line[0] );
    mols.push_back( make_pair( split_line[0] , mol_name ) );
    ++in; // get past '\n'
    ++mol_num;
  }

  cout << "Read " << mols.size() << " SMILES strings." << endl;

}

// ****************************************************************************
void MolDiffViewer::find_differences() {

  for( int i = 0 , is = std::min( file1_mols_.size() , file2_mols_.size() ) ; i < is ; ++i ) {
#ifdef NOTYET
    cout << file1_mols_[i].second << " : " << file2_mols_[i].second << " :: "
         << file1_mols_[i].first << " : " << file2_mols_[i].first << endl;
#endif
    if( file1_mols_[i].second == file2_mols_[i].second ) {
      if( file1_mols_[i].first != file2_mols_[i].first ) {
        diffs_.push_back( i );
      }
    } else {
      break; // because that's as far as we can go
    }
  }

  cout << "Number of differences : " << diffs_.size() << endl;

}

// ****************************************************************************
void MolDiffViewer::show_other_smiles( const string &mol_name ,
                                       const vector<pair<string, string> > &mols ,
                                       QTextEdit *other_smiles ) {

  using namespace boost;
  other_smiles->clear();
  vector<pair<string,string> >::const_iterator p = std::find_if( mols.begin() , mols.end() ,
                                                                 bind( std::equal_to<string>() ,
                                                                       bind( &pair<string,string>::second , _1 ) ,
                                                                       mol_name ) );
  while( p != mols.end() && p->second == mol_name ) {
    other_smiles->insertPlainText( QString( "%1 %2\n" ).arg( p->first.c_str() ).arg( p->second.c_str() ) );
    ++p;
  }

}

// ****************************************************************************
void MolDiffViewer::slot_quit() {

  exit( 0 );

}

// ****************************************************************************
void MolDiffViewer::slot_slider_changed() {

  int diff = diffs_[mol_slider_->value()];

  OEMolBase *mol_l = OENewMolBase( OEMolBaseType::OEDefault );
  OEParseSmiles( *mol_l , file1_mols_[diff].first );
  mol_l->SetTitle( file1_mols_[diff].second );
  mol_disp_l_->set_display_molecule( mol_l );
  delete mol_l;
  smiles_1_->setText( file1_mols_[diff].first.c_str() );

  OEMolBase *mol_r = OENewMolBase( OEMolBaseType::OEDefault );
  OEParseSmiles( *mol_r , file2_mols_[diff].first );
  mol_r->SetTitle( file2_mols_[diff].second );
  mol_disp_r_->set_display_molecule( mol_r );
  delete mol_r;
  smiles_2_->setText( file2_mols_[diff].first.c_str() );

  if( !file3_mols_.empty() || !file4_mols_.empty() ) {
    show_other_smiles( file1_mols_[diff].second , file3_mols_ , other_smiles_1_ );
    show_other_smiles( file1_mols_[diff].second , file4_mols_ , other_smiles_2_ );
  }

  statusBar()->showMessage( QString( "Diff %1 of %2" ).arg( mol_slider_->value() + 1 ).arg( mol_slider_->maximum() + 1 ) );
}
