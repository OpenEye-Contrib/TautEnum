//
// file QTMolDisplay2D.cc
// Dave Cosgrove
// AstraZeneca
// 1st March 2007
//
// This class uses Qt and OEDepict to display an OEMolBase object in a QPainter.

#include <QAction>
#include <QApplication>
#include <QContextMenuEvent>
#include <QMenu>
#include <QPaintEvent>
#include <QToolTip>

#include "stddefs.H"
#include "DACOEMolAtomIndex.H"
#include "DACOEMolBondIndex.H"
#include "QTMolDisplay2D.H"

#include <oedepict.h>

#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OEDepict;
using namespace OESystem;

namespace DACLIB {
// in eponymous file
void split_smiles_into_atom_bits( const string &smiles ,
                                  vector<string> &bits );
string extract_smarts_from_smirks( const std::string &smirks_string );
// in eponymous file
pair<QImage *,OE2DMolDisplay *> draw_oemol_to_qimage( QWidget *wid , OEMolBase &mol ,
                                                      bool coloured_mol ,
                                                      const vector<pair<OEAtomBase * , string> > &atom_labels ,
                                                      const vector<pair<OEAtomBase * , QColor> > &atom_colours ,
                                                      const vector<pair<OEBondBase * , QColor> > &bond_colours );
}

namespace DACLIB {

  // *******************************************************************************
  QTMolDisplay2D::QTMolDisplay2D( QWidget *p , Qt::WindowFlags f ) :
    QWidget( p , f ) , disp_mol_( 0 ) , coloured_mol_( true ) ,
    min_font_size_( 6 ) , line_width_( 1 ) , background_colour_( QColor( "White" ) ) {

    build_actions();

    setMouseTracking( true ); // for the tooltips

  }

  // *******************************************************************************
  QTMolDisplay2D::~QTMolDisplay2D() {

    delete disp_mol_;

  }

  // *******************************************************************************
  void QTMolDisplay2D::set_display_molecule( OEMolBase *new_mol ) {

    if( disp_mol_ ) {
      delete disp_mol_;
    }
    atom_labels_.clear();
    atom_tooltips_.clear();
    atom_colours_.clear();
    sel_atoms_.clear();

    if( new_mol ) {
      disp_mol_ = OENewMolBase( *new_mol , OEMolBaseType::OEDefault );
    } else {
      disp_mol_ = 0;
    }

    if( disp_mol_ ) {
      OEPrepareDepiction( *disp_mol_ , true );
    }

    if( toggle_atom_nums_->isChecked() ) {
      number_atoms();
    }
    update();

  }

  // *******************************************************************************
  void QTMolDisplay2D::clear_display_molecule() {

    if( disp_mol_ ) {
      disp_mol_->Clear();
      atom_labels_.clear();
      atom_tooltips_.clear();
      atom_colours_.clear();
      sel_atoms_.clear();
      update();
    }

  }

  // *******************************************************************************
  // put sequence numbers on the atoms
  void QTMolDisplay2D::number_atoms() {

    atom_labels_.clear();
    OEIter<OEAtomBase> atom;
    for( atom = disp_mol_->GetAtoms() ; atom ; ++atom ) {
    
      string atnam;
      if( !atom->GetAtomicNum() ) {
        if( atom->GetMapIdx() ) {
          atnam = string( "R:" ) += lexical_cast<string>( atom->GetMapIdx() );
        } else {
          atnam = string( "*" );
        }
      }
      atnam += lexical_cast<string>( DACLIB::atom_index( *atom ) + 1 );
      atom_labels_.push_back( make_pair( atom , atnam ) );
    }

  }

  // *******************************************************************************
  // put sequence numbers on the atoms with atomic number 0 and no map index
  // (which will have been input as * in the SMILES string, probably)
  void QTMolDisplay2D::number_star_atoms() {

    atom_labels_.clear();
    int i = 1;
    OEIter<OEAtomBase> atom;
    for( atom = disp_mol_->GetAtoms() ; atom ; ++atom , ++i ) {
      if( !atom->GetAtomicNum() && !atom->GetMapIdx() ) {
        string atnam = lexical_cast<string>( i );
        atom_labels_.push_back( make_pair( atom , atnam ) );
      }
    }

    update();

  }

  // *******************************************************************************
  void QTMolDisplay2D::number_atoms_by_default( bool na ) {

    toggle_atom_nums_->setChecked( na );

  }

  // *******************************************************************************
  void QTMolDisplay2D::colour_atoms( const vector<OEAtomBase *> &atoms ,
                                     const QColor &colour ) {

    // make sure that all the requested atoms are in the molecule. Mayhem can ensue
    // otherwise.
    atom_colours_.clear();
    for( size_t i = 0 , is = atoms.size() ; i < is ; ++i ) {
      OEIter<OEAtomBase> ma =
          disp_mol_->GetAtoms( HasAtomIndex( atom_index( *atoms[i] ) ) );
      if( ma ) {
        atom_colours_.push_back( make_pair( ma , colour ) );
      }
    }

    update();

  }

  // *******************************************************************************
  // colour the atoms whose DAC_atom_index values are given.
  void QTMolDisplay2D::colour_atoms( const vector<unsigned int> &atom_idxs ,
                                     const QColor &colour ) {

    atom_colours_.clear();
    for( size_t i = 0 , is = atom_idxs.size() ; i < is ; ++i ) {
      OEIter<OEAtomBase> ats = disp_mol_->GetAtoms( HasAtomIndex( atom_idxs[i] ) );
      if( ats ) {
        atom_colours_.push_back( make_pair( ats , colour ) );
      }
    }

    update();

  }

  // *******************************************************************************
  void QTMolDisplay2D::colour_atoms( const QColor &colour ) {

    OEIter<OEAtomBase> atom;
    vector<OEAtomBase *> atoms;
    atoms.reserve( disp_mol_->NumAtoms() );
    for( atom = disp_mol_->GetAtoms() ; atom ; ++atom )
      atoms.push_back( atom );

    colour_atoms( atoms , colour );

  }

  // *******************************************************************************
  // returns the number of atoms in the SMARTS hit
  int QTMolDisplay2D::colour_atoms( const string &smarts_string ,
                                    const QColor &colour ) {

    OESubSearch subsearch( smarts_string.c_str() );
    OEIter<OEMatchBase> match = subsearch.Match( *disp_mol_ , true );

    vector<OEAtomBase *> match_ats;
    int num_atoms_in_smarts_hit = 0;

    for( ; match ; ++match ) {
      OEIter<OEAtomBase> ma;
      num_atoms_in_smarts_hit = 0;
      for( ma = match->GetTargetAtoms() ; ma ; ++ma ) {
        match_ats.push_back( ma );
        ++num_atoms_in_smarts_hit;
      }
    }

    colour_atoms( match_ats , colour );

    return num_atoms_in_smarts_hit;

  }

  // *******************************************************************************
  // colour atoms hit by the OESubSearch objects. Atoms hit once will be Red,
  // twice Orange, then Yellow, Green, Blue and Purple for > 5.
  void QTMolDisplay2D::colour_atoms( const vector<pair<boost::shared_ptr<OESubSearch>,string> > &sub_searches ) {

    vector<int> hit_counts( disp_mol_->GetMaxAtomIdx() , 0 );

    atom_tooltips_.clear();

    for( size_t i = 0 , is = sub_searches.size() ; i < is ; ++i ) {
      // unique matches only (that's the true bit in Match)
      OEIter<OEMatchBase> match = sub_searches[i].first->Match( *disp_mol_ , true );
      for( ; match ; ++match ) {
        for( OEIter<OEAtomBase> ma = match->GetTargetAtoms() ; ma ; ++ma ) {
          append_atom_tooltip( ma , sub_searches[i].second );
          ++hit_counts[ma->GetIdx()];
        }
      }
    }

    static vector<QColor> count_colours;
    if( count_colours.empty() ) {
      count_colours.push_back( QColor( "Black") ); // shouldn't ever be used
      count_colours.push_back( QColor( "Red") );
      count_colours.push_back( QColor( "Orange") );
      count_colours.push_back( QColor( "Yellow") );
      count_colours.push_back( QColor( "Green") );
      count_colours.push_back( QColor( "Blue") );
      count_colours.push_back( QColor( "Purple") );
    }

    atom_colours_.clear();
    for( OEIter<OEAtomBase> at = disp_mol_->GetAtoms() ; at ; ++at ) {
      unsigned int at_ind = at->GetIdx();
      if( hit_counts[at_ind] > 0 && hit_counts[at_ind] < 6 ) {
        atom_colours_.push_back( make_pair( at , count_colours[hit_counts[at_ind]] ) );
      } else if( hit_counts[at_ind] > 5 ) {
        atom_colours_.push_back( make_pair( at , count_colours[5] ) );
      }
    }

    update();

  }

  // *******************************************************************************
  // colour the bonds passed in, using the DACLIB::bond_index() values.
  void QTMolDisplay2D::colour_bonds( const std::vector<unsigned int> &bond_idxs ,
                                     const QColor &colour ) {

    bond_colours_.clear();
    for( size_t i = 0 , is = bond_idxs.size() ; i < is ; ++i ) {
      OEIter<OEBondBase> bds = disp_mol_->GetBonds( DACLIB::HasBondIndex( bond_idxs[i] ) );
      if( bds ) {
        bond_colours_.push_back( make_pair( bds , colour ) );
      }
    }

  }

  // *******************************************************************************
  // put sequence numbers on the atoms hit by the SMARTS string
  void QTMolDisplay2D::label_atoms_by_smarts( const std::string &smarts_string ) {

    OESubSearch subsearch( smarts_string.c_str() );
    OEIter<OEMatchBase> match = subsearch.Match( *disp_mol_ , true );

    atom_labels_.clear();
    for( ; match ; ++match ) {
      OEIter<OEAtomBase> mp;
      int i = 1;
      for( mp = match->GetTargetAtoms() ; mp ; ++mp , ++i ) {
        string atnam;
        if( mp->GetAtomicNum() )
          atnam = string( OEGetAtomicSymbol( mp->GetAtomicNum() ) );
        else {
          if( mp->GetMapIdx() )
            atnam = string( "R:" ) + lexical_cast<string>( mp->GetMapIdx() );
          else
            atnam = string( "*" );
        }
        atnam += lexical_cast<string>( i );
        atom_labels_.push_back( make_pair( mp , atnam ) );
      }
    }
  
    update();

  }

// *******************************************************************************
// label the atoms by sequence number of the atoms in vector order
void QTMolDisplay2D::label_atoms_by_number( const std::vector<OEChem::OEAtomBase *> &atoms ) {

  atom_labels_.clear();
  atom_tooltips_.clear();
  for( size_t i = 0 , is = atoms.size() ; i < is ; ++i ) {
    OEIter<OEAtomBase> ma =
      disp_mol_->GetAtoms( DACLIB::HasAtomIndex( DACLIB::atom_index( *atoms[i] ) ) );
    string atnam( "H" );
    if( ma->GetAtomicNum() ) {
      atnam = lexical_cast<string>( i + 1 );
    }
    atom_labels_.push_back( make_pair( ma , atnam ) );
  }

  update();

}

// *******************************************************************************
// label the atoms by a MapIdx if they have them
void QTMolDisplay2D::label_atoms_by_map_idx() {

  atom_labels_.clear();
  atom_tooltips_.clear();
  for( OEIter<OEAtomBase> atom = disp_mol_->GetAtoms() ; atom ; ++atom ) {
    string atnam;
    if( atom->GetMapIdx() ) {
      if( atom->GetAtomicNum() ) {
        atnam = OEGetAtomicSymbol( atom->GetAtomicNum() );
      }
      atnam += lexical_cast<string>( atom->GetMapIdx() );
      atom_labels_.push_back( make_pair( atom , atnam ) );
    }
  }

  update();

}

// *******************************************************************************
// this one uses the first half of the SMIRKS string to supply the
// the atom map indices by parsing them out directly
// returns the SMARTS string extracted and used.
void QTMolDisplay2D::label_atoms_by_map_idx( const string &smirks_string ,
                                             string &smarts_string ) {

  smarts_string = DACLIB::extract_smarts_from_smirks( smirks_string );

  vector<string> atom_bits;
  split_smiles_into_atom_bits( smarts_string , atom_bits );
#ifdef NOTYET
  cout << "SMARTS : " << smarts_string << endl;
  cout << "atom bits : ";
  copy( atom_bits.begin() , atom_bits.end() , stringOut );
  cout << endl;
#endif

  // extract the map_idxs from the atom_bits
  vector<unsigned int> map_idxs( atom_bits.size() , 0 );
  for( size_t k = 0 , ks = atom_bits.size() ; k < ks ; ++k ) {
    size_t i = atom_bits[k].rfind( ":" );
    if( i == string::npos ) {
      continue; // not index for this atom
    }
    size_t j = atom_bits[k].find( ']' , i );
    string sidx( atom_bits[k].substr( i + 1 , j - i - 1 ) );
    try {
      unsigned int idx = lexical_cast<unsigned int>( sidx ) ;
      map_idxs[k] = idx;
    } catch( bad_lexical_cast &e ) {
      // don't need to do anything except catch the exception
      // if, for some reason, there's not a number between the
      // last : and the end of the atom_bit
    }
  }

  OESubSearch subsearch( smarts_string.c_str() );
  OEIter<OEMatchBase> match = subsearch.Match( *disp_mol_ , true );

  // set the MapIdx properties for all the relevant atoms
  for( ; match ; ++match ) {
    OEIter<OEAtomBase> mp;
    int i = 0;
    for( mp = match->GetTargetAtoms() ; mp ; ++mp , ++i ) {
      if( map_idxs[i] ) {
        mp->SetMapIdx( map_idxs[i] );
      }
    }
  }

  label_atoms_by_map_idx();

}

// *******************************************************************************
// select the atoms according to the SMARTS string passed in
void QTMolDisplay2D::select_atoms_by_smarts( const std::string &smarts_string ) {

  cout << "QTMolDisplay2D::select_atoms_by_smarts : " << smarts_string << endl;

  OESubSearch subsearch( smarts_string.c_str() );
  OEIter<OEMatchBase> match = subsearch.Match( *disp_mol_ , true );

  sel_atoms_.clear();
  for( ; match ; ++match ) {
    OEIter<OEAtomBase> mp;
    for( mp = match->GetTargetAtoms() ; mp ; ++mp ) {
      sel_atoms_.push_back( mp );
    }
  }

  update();

}

// *******************************************************************************
void QTMolDisplay2D::set_line_width( int new_width ) {

  line_width_ = new_width;
  update();

}

// *******************************************************************************
  void QTMolDisplay2D::set_background_colour( const QColor &colour ) {

    background_colour_ = colour;
    update();

  }

  // *******************************************************************************
  void QTMolDisplay2D::build_actions() {

    toggle_atom_nums_ = new QAction( "Show Numbers" , this );
    toggle_atom_nums_->setCheckable( true );
    toggle_atom_nums_->setChecked( false );
    connect( toggle_atom_nums_ , SIGNAL( triggered() ) ,
             this , SLOT( slot_toggle_atom_nums() ) );

    toggle_black_white_ = new QAction( "Draw Black+White" , this );
    toggle_black_white_->setCheckable( true );
    toggle_black_white_->setChecked( false );
    connect( toggle_black_white_ , SIGNAL( triggered() ) ,
             this , SLOT( slot_toggle_black_white() ) );

  }

  // *******************************************************************************
  void QTMolDisplay2D::render_atom_labels() {

    for( size_t i = 0 , is = atom_labels_.size() ; i < is ; ++i ) {
      OE2DAtomDisplay *adisp = disp_->GetAtomDisplay( atom_labels_[i].first );
      adisp->SetProperty( atom_labels_[i].second );
    }

  }

  // *************************************************************************
  void QTMolDisplay2D::squares_round_selected_atoms( QImage *image ) {

    QPainter qp( image );
    qp.initFrom( this );
    qp.setRenderHint( QPainter::Antialiasing , true );

    int rect_size = int( width() / 100 );
    if( rect_size < 5 ) {
      rect_size = 5;
    }

    for( size_t i = 0 , is = sel_atoms_.size() ; i < is ; ++i ) {
      OE2DAtomDisplay *adisp = disp_->GetAtomDisplay( sel_atoms_[i] );
      OE2DPoint cds = adisp->GetCoords();
      qp.setPen( "Orange" );
      qp.drawRect( int( cds.GetX() ) - rect_size , int( cds.GetY() ) - rect_size ,
                   2 * rect_size , 2 * rect_size );
    }

  }

  // ************************************************************************
  bool QTMolDisplay2D::event( QEvent *e ) {

    if( e->type() == QEvent::ToolTip ) {
      QHelpEvent *help_event = static_cast<QHelpEvent *>(e);
      show_atom_tooltip( help_event );
    }
    return QWidget::event( e );

  }

  // *******************************************************************************
  void QTMolDisplay2D::paintEvent( QPaintEvent *e __attribute__((unused)) ) {

#ifdef NOTYET
    cout << "QTMolDisplay2D::paintEvent" << endl;
#endif

    QImage *mol_img = draw_molecule();

    QPainter wp( this );
    wp.drawImage( 0 , 0 , *mol_img );

    delete mol_img;

  }

  // *******************************************************************************
  void QTMolDisplay2D::mousePressEvent( QMouseEvent *e ) {

    if( Qt::LeftButton == e->button() ) {
      OEAtomBase *pa = find_nearest_atom( e->x() , e->y() );
      if( !pa && Qt::ControlModifier != e->modifiers() ) {
        sel_atoms_.clear();
      }
      if( pa ) {
        sel_atoms_.push_back( pa );
        emit atom_selected( pa->GetIdx() );
      }
      update();
    }

  }

  // *******************************************************************************
  void QTMolDisplay2D::contextMenuEvent( QContextMenuEvent *e ) {

    QMenu menu( this );
    menu.addAction( toggle_atom_nums_ );
    menu.addAction( toggle_black_white_ );

    menu.exec( e->globalPos() );

  }

  // *******************************************************************************
  QImage *QTMolDisplay2D::draw_molecule() {

    if( !disp_mol_ || !(*disp_mol_) ) {
      QImage *ret_val = new QImage( size() , QImage::Format_ARGB32_Premultiplied );
      QPainter qp( ret_val );
      QColor back_colour( QApplication::palette().color( QPalette::Normal ,
                                                         QPalette::Base ) );

      qp.setBackground( back_colour );
      qp.fillRect( rect() , back_colour );
      return ret_val;
    } else {
      pair<QImage *,OE2DMolDisplay *> img_mol_disp = draw_oemol_to_qimage( this , *disp_mol_ ,
                                                                           coloured_mol_ , atom_labels_ ,
                                                                           atom_colours_ , bond_colours_ );

      disp_.reset( img_mol_disp.second );

      // add selected atom squares - must be done after molecule rendering
      squares_round_selected_atoms( img_mol_disp.first );

      return img_mol_disp.first;
    }

  }

  // *************************************************************************
  OEAtomBase *QTMolDisplay2D::find_nearest_atom( int x_pos , int y_pos ) const {

    int nearest_dist = numeric_limits<int>::max();
    OEAtomBase *nearest_atom = 0;

    if( !disp_mol_ || !disp_ ) {
      return nearest_atom;
    }

    int miss_dist = DACLIB::square( int( width() / 50 ) );
    if( miss_dist < 4 ) {
      miss_dist = 4;
    }

    OEIter<OEAtomBase> atom;
    for( atom = disp_mol_->GetAtoms() ; atom ; ++atom ) {
      OE2DAtomDisplay *adisp = disp_->GetAtomDisplay( atom );
      OE2DPoint cds = adisp->GetCoords();

      int sq_dist = DACLIB::square( int( cds.GetX() ) - x_pos ) + DACLIB::square( int( cds.GetY() ) - y_pos );
      if( sq_dist < nearest_dist && sq_dist < miss_dist ) {
        nearest_atom = atom;
        nearest_dist = sq_dist;
      }
    }

    return nearest_atom;

  }

  // *******************************************************************************
  // append the given string to the tooltip of the atom
  void QTMolDisplay2D::append_atom_tooltip( OEAtomBase *atom , const string &label ) {

    int i = -1;
    for( int j = 0 , js = static_cast<int>( atom_tooltips_.size() ) ; j < js ; ++j ) {
      if( atom_tooltips_[j].first == atom ) {
        i = j;
        break;
      }
    }
    if( i == -1 ) {
      atom_tooltips_.push_back( make_pair( atom , label.c_str() ) );
    } else {
      atom_tooltips_[i].second += QString( ",%1" ).arg( label.c_str() );
    }

  }

  // *******************************************************************************
  void QTMolDisplay2D::show_atom_tooltip( QHelpEvent *e ) {

    OEAtomBase *atom = find_nearest_atom( e->pos().x() , e->pos().y() );
    if( atom ) {
      for( size_t i = 0 , is = atom_tooltips_.size() ; i < is ; ++i ) {
        if( atom_tooltips_[i].first == atom ) {
          QToolTip::showText( e->globalPos() , atom_tooltips_[i].second );
          break;
        }
      }
    }
  }

  // *******************************************************************************
  void QTMolDisplay2D::slot_atom_selected( unsigned int oe_ind ) {

    if( !disp_mol_ ) {
      return;
    }

    sel_atoms_.clear();
    OEIter<OEAtomBase> at = disp_mol_->GetAtoms( OEHasAtomIdx( oe_ind ) );
    if( at ) {
      sel_atoms_.push_back( at );
      update();
    }
  }

  // *******************************************************************************
  void QTMolDisplay2D::slot_toggle_atom_nums() {

    if( toggle_atom_nums_->isChecked() ) {
      number_atoms();
      update();
    } else {
      clear_atom_labels();
      update();
    }

  }

  // *******************************************************************************
  void QTMolDisplay2D::slot_toggle_black_white() {

    set_coloured_mol( !toggle_black_white_->isChecked() );
    update();

  }

} // end of namespace DACLIB
