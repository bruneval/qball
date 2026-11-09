////////////////////////////////////////////////////////////////////////////////
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// DensityDipole.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef DENSITYDIPOLE_H
#define DENSITYDIPOLE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

// Print the electronic dipole computed from the density,
//   <electronic_dipole_density> -int r n(r) d3r </electronic_dipole_density>
// at every RT-TDDFT step, with r taken on the cell grid (r = i/np0 a + ...).
// Meaningful only for an isolated system whose density is negligible at the
// cell boundaries. Default OFF.

class DensityDipole : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "density_dipole"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " density_dipole takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " density_dipole must be ON or OFF" << endl;
      return 1;
    }

    s->ctrl.density_dipole = ( v == "ON" );

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) ;
     if ( s->ctrl.density_dipole ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  DensityDipole(Sample *sample) : s(sample) { s->ctrl.density_dipole = false; }
};
#endif

// Local Variables:
// mode: c++
// End:
