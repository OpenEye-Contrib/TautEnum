//
// file QTMolDisplay2D.H
// David Cosgrove
// AstraZeneca
// 1st March 2007
//
// This class uses Qt and OEDepict to display an OEMolBase object in a QPainter.

#ifndef DAC_QT_MOL_DISPLAY_2D
#define DAC_QT_MOL_DISPLAY_2D

#include <vector>

#include <QPainter>
#include <QSize>
#include <QString>
#include <QWidget>

#include <oechem.h>
#include <oedepict.h>

#include <boost/shared_ptr.hpp>

// ******************************************************************************

class QAction;
class QColor;
class QHelpEvent;
class QMouseEvent;
class QPainter;
namespace OEChem {
  class OESubSearch;
}

// ******************************************************************************

namespace DACLIB {

  class QTMolDisplay2D : public QWidget {

    Q_OBJECT

    signals :

    void atom_selected( unsigned int oe_ind );

    public slots :

    void slot_atom_selected( unsigned int oe_ind );
    void slot_toggle_atom_nums();
    void slot_toggle_black_white();

  protected :

    OEChem::OEMolBase *disp_mol_;
    boost::shared_ptr<OEDepict::OE2DMolDisplay> disp_;

    // atom_labels_ are small bits of text displayed as part of the atom name, such
    // as the sequence number from a SMARTS pattern
    std::vector<std::pair<OEChem::OEAtomBase * , std::string> > atom_labels_;
    // atom_tool_tips_ are longer bits of text shown as a tooltip when the user
    // hovers over the atom
    std::vector<std::pair<OEChem::OEAtomBase * , QString> > atom_tooltips_;
    std::vector<std::pair<OEChem::OEAtomBase * , QColor> > atom_colours_;
    std::vector<std::pair<OEChem::OEBondBase * , QColor> > bond_colours_;
    std::vector<OEChem::OEAtomBase *> sel_atoms_;

    bool coloured_mol_; /* black on white or coloured atom symbols - not to be confused with output
                           of colour_atoms() and colour_bonds() routines, which take precedence */
    int min_font_size_; // -1 (default) means no minimum
    int line_width_; // for the lines in the molecule display

  private :

    QAction *toggle_atom_nums_ , *toggle_black_white_;
    QColor background_colour_;

    void build_actions();

    void render_atom_labels();
    void squares_round_selected_atoms( QImage *image );

  protected :

    bool event( QEvent *e ); // overloaded for tooltips
    void paintEvent( QPaintEvent *e );
    void mousePressEvent( QMouseEvent *e );
    void contextMenuEvent( QContextMenuEvent *e );

    virtual QImage *draw_molecule();
    // nearest atom to screen coords
    OEChem::OEAtomBase *find_nearest_atom( int x_pos , int y_pos ) const;

    // append the given string to the tooltip of the atom
    void append_atom_tooltip( OEChem::OEAtomBase *atom , const std::string &label );
    void show_atom_tooltip( QHelpEvent *e );

  public :

    QTMolDisplay2D( QWidget *p = 0 , Qt::WindowFlags f = 0 );
    virtual ~QTMolDisplay2D();

    QSize sizeHint() const { return QSize( 200 , 200 ); }

    std::vector<OEChem::OEAtomBase *> sel_atoms() const {
      return sel_atoms_;
    }

    // takes a copy of new_mol, so the one passed in doesn't need to last
    virtual void set_display_molecule( OEChem::OEMolBase *new_mol );
    virtual void clear_display_molecule();
    OEChem::OEMolBase *display_molecule() { return disp_mol_; }
    boost::shared_ptr<OEDepict::OE2DMolDisplay> oemoldisplay() { return disp_; }

    // put sequence numbers on the atoms
    void number_atoms();
    // put sequence numbers on the atoms with atomic number 0 and no map index
    // (which will have been input as * in the SMILES string, probably)
    void number_star_atoms();
    void number_atoms_by_default( bool na = true );

    void colour_atoms( const std::vector<OEChem::OEAtomBase *> &atoms ,
                       const QColor &colour );
    // colour the atoms whose DACLIB::atom_index() values are given.
    void colour_atoms( const std::vector<unsigned int> &atom_idxs ,
                       const QColor &colour );
    void colour_atoms( const QColor &colour );
    // next overload returns the number of atoms in the SMARTS hit
    int colour_atoms( const std::string &smarts_string , const QColor &colour );
    // colour atoms hit by the OESubSearch objects. Atoms hit once will be Red,
    // twice Orange, then Yellow, Green, Blue and Purple for > 5.
    void colour_atoms( const std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > &sub_searches );

    // colour the bonds passed in, using the DACLIB::bond_index() values.
    void colour_bonds( const std::vector<unsigned int> &bond_idx ,
                       const QColor &colour );

    // put sequence numbers on the atoms hit by the SMARTS string
    void label_atoms_by_smarts( const std::string &smarts_string );
    // label the atoms by sequence number of the atoms in vector order
    void label_atoms_by_number( const std::vector<OEChem::OEAtomBase *> &atoms );
    // label the atoms by a MapIdx if they have them
    void label_atoms_by_map_idx();
    // this one uses the first half of the SMIRKS string to supply the
    // the atom map indices by parsing them out directly
    // returns the SMARTS string extracted and used.
    void label_atoms_by_map_idx( const std::string &smirks_string ,
                                 std::string &smarts_string );

    // select the atoms according to the SMARTS string passed in
    void select_atoms_by_smarts( const std::string &smarts_string );

    void clear_atom_labels() {
      atom_labels_.clear();
      update();
    }
    void clear_atom_colours() {
      atom_colours_.clear();
      update();
    }
    void set_coloured_mol( bool new_val ) {
      if( new_val != coloured_mol_ ) {
        coloured_mol_ = new_val;
        update();
      }
    }
    int min_font_size() const { return min_font_size_; }
    void set_min_font_size( int new_val ) { min_font_size_ = new_val; }
    void set_line_width( int new_width );

    void set_background_colour( const QColor &colour );

  };

}

#endif

