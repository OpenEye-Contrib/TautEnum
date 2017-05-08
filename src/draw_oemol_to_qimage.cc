//
// file draw_mol_to_qimage.cc
// David Cosgrove
// AstraZeneca
// 2nd April 2014
//
// Takes an OEMol and various other things and returns a QImage with the
// the molecule drawn in it.

#include <QApplication>
#include <QImage>
#include <QPainter>
#include <QString>
#include <QWidget>

#include <oechem.h>
#include <oedepict.h>

#include <vector>

#include "DACOEMolAtomIndex.H"

using namespace std;
using namespace OEChem;
using namespace OEDepict;
using namespace OESystem;

namespace DACLIB {

// **********************************************************************
void render_atom_labels( OE2DMolDisplay &mol_disp ,
                         const vector<pair<OEAtomBase * , string> > &atom_labels ) {

  for( size_t i = 0 , is = atom_labels.size() ; i < is ; ++i ) {
    mol_disp.GetAtomDisplay( atom_labels[i].first )->SetProperty( atom_labels[i].second );
  }

}

// *************************************************************************
void add_atom_colours_to_depiction( OE2DMolDisplay &mol_disp ,
                                    const vector<pair<OEAtomBase * , QColor> > &atom_colours ) {

  for( size_t i = 0 , is = atom_colours.size() ; i < is ; ++i ) {
    OE2DAtomDisplay *adisp = mol_disp.GetAtomDisplay( atom_colours[i].first );
    OEFont font = adisp->GetLabelFont();
    font.SetColor( OEColor( atom_colours[i].second.red() ,
			    atom_colours[i].second.green() ,
			    atom_colours[i].second.blue() ) );
    adisp->SetLabelFont( font );

    for( OEIter<OEBondBase> bond = atom_colours[i].first->GetBonds() ; bond ; ++bond ) {
      OE2DBondDisplay *bdisp = mol_disp.GetBondDisplay( bond );
      if( bond->GetBgn() == atom_colours[i].first ) {
        OEPen pen = bdisp->GetBgnPen();
        pen.SetForeColor( OEColor( atom_colours[i].second.red() ,
				   atom_colours[i].second.green() ,
				   atom_colours[i].second.blue() ) );
        bdisp->SetBgnPen( pen );
      } else {
        OEPen pen = bdisp->GetEndPen();
        pen.SetForeColor( OEColor( atom_colours[i].second.red() ,
				   atom_colours[i].second.green() ,
				   atom_colours[i].second.blue() ) );
        bdisp->SetEndPen( pen );
      }
    }
  }

}

// *************************************************************************
void add_bond_colours_to_depiction( OE2DMolDisplay &mol_disp ,
                                    const vector<pair<OEBondBase * , QColor> > &bond_colours ) {

  for( size_t i = 0 , is = bond_colours.size() ; i < is ; ++i ) {
    OE2DBondDisplay *bdisp = mol_disp.GetBondDisplay( bond_colours[i].first );
    OEPen pen = bdisp->GetBgnPen();
    pen.SetForeColor( OEColor( bond_colours[i].second.red() ,
			       bond_colours[i].second.green() ,
			       bond_colours[i].second.blue() ) );
    pen = bdisp->GetEndPen();
    pen.SetForeColor( OEColor( bond_colours[i].second.red() ,
			       bond_colours[i].second.green() ,
			       bond_colours[i].second.blue() ) );
  }

}

// ****************************************************************************
OE2DMolDisplay *create_oe_mol_display( QWidget *wid , OEMolBase &mol ,
                                       bool coloured_mol ,
                                       const vector<pair<OEAtomBase * , string> > &atom_labels ,
                                       const vector<pair<OEAtomBase * , QColor> > &atom_colours ,
                                       const vector<pair<OEBondBase * , QColor> > &bond_colours ) {

  OEPrepareDepictionOptions popts;
  OEPrepareDepiction( mol , popts );

  OE2DMolDisplayOptions dopts( double( wid->width() ) , double( wid->height() ) ,
                               OEScale::AutoScale );
  dopts.SetTitleLocation( OETitleLocation::Bottom );
  dopts.SetBondStereoStyle( OEBondStereoStyle::Display::CIPBondStereo );
  OEPen pen;
  pen.SetLineWidth( 1.0 );
  dopts.SetDefaultBondPen( pen );
  if( !coloured_mol ) {
    dopts.SetAtomColorStyle( OEAtomColorStyle::WhiteMonochrome );
  } else {
    dopts.SetAtomColorStyle( OEAtomColorStyle::WhiteCPK );
  }

  dopts.SetTitleFontScale( 0.5 );
  OE2DMolDisplay *mol_disp = new OE2DMolDisplay( mol , dopts );

  if( !coloured_mol ) {
    add_atom_colours_to_depiction( *mol_disp , atom_colours );
    add_bond_colours_to_depiction( *mol_disp , bond_colours );
  }

  // add the atom labels
  render_atom_labels( *mol_disp , atom_labels );

  return mol_disp;

}

// ****************************************************************************
QImage *oe_mol_disp_to_qimage( OE2DMolDisplay &mol_disp ) {

  OEPlatform::oeosstream oes;
  OERenderMolecule( oes , "PNG" , mol_disp );

  QImage *img = new QImage;
  img->loadFromData( reinterpret_cast<const unsigned char *>( oes.str().c_str() ) ,
		     static_cast<int>(oes.str().length()) , "png" );

  return img;

}

// ****************************************************************************
pair<QImage *,OE2DMolDisplay *> draw_oemol_to_qimage( QWidget *wid , OEMolBase &mol ,
                                                      bool coloured_mol ,
                                                      const vector<pair<OEAtomBase * , string> > &atom_labels ,
                                                      const vector<pair<OEAtomBase * , QColor> > &atom_colours ,
                                                      const vector<pair<OEBondBase * , QColor> > &bond_colours ) {


  OE2DMolDisplay *mol_disp = create_oe_mol_display( wid , mol , coloured_mol ,
                                                    atom_labels , atom_colours ,
                                                    bond_colours );
  QImage *img = oe_mol_disp_to_qimage( *mol_disp );

  return make_pair( img , mol_disp );

}

} // EO namespace DACLIB
