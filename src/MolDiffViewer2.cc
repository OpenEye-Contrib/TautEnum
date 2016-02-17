//
// file MolDiffViewer2.cc
// David Cosgrove
// AstraZeneca
// 22nd November 2013

#include "MolDiffViewer2.H"

#include "QTMolDisplay2D.H"

#include <QAction>
#include <QFileDialog>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QSlider>
#include <QSplitter>
#include <QStatusBar>
#include <QString>
#include <QTextEdit>

#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>

#include <oechem.h>

using namespace OEChem;
using namespace boost;
using namespace std;

// ****************************************************************************
MolDiffViewer2::MolDiffViewer2( int argc , char **argv ) : QMainWindow() {

  if( argc < 2 ) {
    cerr << "MolDiffViewer2 : needs the name of 2 SMILES files." << endl;
    exit( 1 );
  }
  build_widget();

  read_smiles_file( argv[1] , file1_mols_ );
  file1_label_->setText( argv[1] );
  read_smiles_file( argv[2] , file2_mols_ );
  file2_label_->setText( argv[2] );

  find_differences();
  if( !diffs_.empty() ) {
    mol_slider_->setMaximum( diffs_.size() - 1 );
    mol_slider_->setMinimum( 0 );
    slot_slider_changed();
  } else {
    mol_slider_->setEnabled( false );
  }

}

// ****************************************************************************
void MolDiffViewer2::build_widget() {

  build_actions();
  build_menubar();

  mol_slider_ = new QSlider;
  mol_slider_->setPageStep( 1 );
  connect( mol_slider_ , SIGNAL( valueChanged(int) ) ,
           this , SLOT( slot_slider_changed() ) );

  QVBoxLayout *vbox1 = new QVBoxLayout;
  QSplitter *s1 = new QSplitter( Qt::Vertical );
  QWidget *w1 = new QWidget;
  left_grid_ = new QGridLayout;
  w1->setLayout( left_grid_ );
  s1->addWidget( w1 );
  smiles_1_ = new QTextEdit;
  s1->addWidget( smiles_1_ );
  vbox1->addWidget( s1 );
  file1_label_ = new QLabel;
  vbox1->addWidget( file1_label_ );

  QVBoxLayout *vbox2 = new QVBoxLayout;
  QSplitter *s2 = new QSplitter( Qt::Vertical );
  QWidget *w2 = new QWidget;
  right_grid_ = new QGridLayout;
  w2->setLayout( right_grid_ );
  s2->addWidget( w2 );
  smiles_2_ = new QTextEdit;
  s2->addWidget( smiles_2_ );
  vbox2->addWidget( s2 );
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
void MolDiffViewer2::build_actions() {

  file_quit_ = new QAction( "Quit" , this );
  file_quit_->setShortcut( QString( "Ctrl+Q" ) );
  connect( file_quit_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_quit() ) );

  next_left_smaller_ = new QAction( "Next Left Smaller" , this );
  next_left_smaller_->setShortcut( QString( "Ctrl+S" ) );
  connect( next_left_smaller_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_next_left_smaller() ) );
  next_left_smaller_->setToolTip( QString( "Goes to next difference where number of molecules smaller in left panel than right." ) );

  write_diffs_ = new QAction( "Write Diffs" , this );
  connect( write_diffs_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_write_diffs() ) );
  write_diffs_->setToolTip( QString( "Write names of molecules that have differences in two files." ) );

}

// ****************************************************************************
void MolDiffViewer2::build_menubar() {

  QMenu *file_menu = menuBar()->addMenu( "File" );
  file_menu->addAction( next_left_smaller_ );
  file_menu->addAction( write_diffs_ );
  file_menu->addAction( file_quit_ );

}

// ****************************************************************************
void MolDiffViewer2::read_smiles_file( const string &filename ,
                                       vector<pair<string,vector<string> > > &mols ) {

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
    if( mols.empty() || mol_name != mols.back().first ) {
      if( !mols.empty() ) {
        sort( mols.back().second.begin() , mols.back().second.end() );
      }
      mols.push_back( make_pair( mol_name , vector<string>( 1 , split_line[0] ) ) );
    } else {
      mols.back().second.push_back( split_line[0] );
    }
    ++in; // get past '\n'
    ++mol_num;
  }

  cout << "Read " << mols.size() << " SMILES strings." << endl;

}

// ****************************************************************************
void MolDiffViewer2::find_differences() {

  for( int i = 0 , is = std::min( file1_mols_.size() , file2_mols_.size() ) ; i < is ; ++i ) {
    if( file1_mols_[i].first != file2_mols_[i].first ) {
      break; // the names need to be in sync
    }
    if( file1_mols_[i].second.size() != file2_mols_[i].second.size() ) {
      diffs_.push_back( i );
    } else {
      // both vectors should be sorted in the same order, so if they're the same
      // they'll have the same elements in the same order
      for( int j = 0 , js = file1_mols_[i].second.size() ; j < js ; ++j ) {
        if( file1_mols_[i].second[j] != file2_mols_[i].second[j] ) {
          diffs_.push_back( i );
          break;
        }
      }
    }
  }

  statusBar()->showMessage( QString( "Number of differences : %1" ).arg( diffs_.size() ) );

}

// ****************************************************************************
void MolDiffViewer2::show_smiles( pair<string, vector<string> > &mols ,
                                  vector<DACLIB::QTMolDisplay2D *> &mol_disps ,
                                  QGridLayout *grid , QTextEdit *smiles ) {

  smiles->clear();
  for( int i = 0 , is = mol_disps.size() ; i < is ; ++i ) {
    mol_disps[i]->hide();
  }

  for( int i = 0 , is = mols.second.size() ; i < is ; ++i ) {
    if( i == static_cast<int>( mol_disps.size() ) ) {
      int r = mol_disps.size() / 2;
      int c = mol_disps.size() % 2;
      mol_disps.push_back( new DACLIB::QTMolDisplay2D );
      grid->addWidget( mol_disps.back() , r , c );
      // cout << "New widget " << i << " into row " << r << " column " << c << endl;
    }
    OEMolBase *mol = OENewMolBase( OEMolBaseType::OEDefault );
    OEParseSmiles( *mol , mols.second[i] );
    mol->SetTitle( mols.first );
    mol_disps[i]->set_display_molecule( mol );
    mol_disps[i]->show();
    delete mol;
    smiles->insertPlainText( QString( "%1\n" ).arg( mols.second[i].c_str() ) );
  }

}

// ****************************************************************************
void MolDiffViewer2::slot_quit() {

  exit( 0 );

}

// ****************************************************************************
void MolDiffViewer2::slot_slider_changed() {

  int diff = diffs_[mol_slider_->value()];

  show_smiles( file1_mols_[diff] , mol_disp_ls_ , left_grid_ , smiles_1_ );
  show_smiles( file2_mols_[diff] , mol_disp_rs_ , right_grid_ , smiles_2_ );

  statusBar()->showMessage( QString( "Diff %1 of %2" ).arg( mol_slider_->value() + 1 ).arg( mol_slider_->maximum() + 1 ) );

}

// ****************************************************************************
void MolDiffViewer2::slot_next_left_smaller() {

  int diff = mol_slider_->value() + 1;
  while( diff < mol_slider_->maximum() ) {
#ifdef NOTYET
    cout << diff << " : " << file1_mols_[diffs_[diff]].second.size() << " and "
            << file2_mols_[diffs_[diff]].second.size() << endl;
#endif
    if( file1_mols_[diffs_[diff]].second.size() < file2_mols_[diffs_[diff]].second.size() ) {
      break;
    }
    ++diff;
  }

  mol_slider_->setValue( diff );

  // slot_slider_changed();

}

// ****************************************************************************
void MolDiffViewer2::slot_write_diffs() {

  QString fn = QFileDialog::getSaveFileName( this , "File for diffs" );
  if( fn.isNull() || fn.isEmpty() ) {
    return;
  }

  ofstream ofs( fn.toLocal8Bit().data() );
  BOOST_FOREACH( int d , diffs_ ) {
    ofs << file1_mols_[d].first << " : ";
    copy( file1_mols_[d].second.begin() , file1_mols_[d].second.end() ,
          ostream_iterator<string>( ofs , " " ) );
    ofs << " : ";
    copy( file2_mols_[d].second.begin() , file2_mols_[d].second.end() ,
          ostream_iterator<string>( ofs , " " ) );
    ofs << endl;
  }

}
