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
// TravelingPseudo.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef TRAVELINGPSEUDO_H
#define TRAVELINGPSEUDO_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

// Traveling pseudopotentials (Stengel, Royo, Artacho, PRL 136, 196401 (2026)):
// in RT-TDDFT with atoms_dyn IMPULSIVE, the nonlocal projectors of each moving
// atom are evaluated at k+G-A-v, where v is the atom velocity.
// Default OFF: results are then identical to the standard implementation.

class TravelingPseudo : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "traveling_pseudo"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " traveling_pseudo takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " traveling_pseudo must be ON or OFF" << endl;
      return 1;
    }

    s->ctrl.traveling_pseudo = ( v == "ON" );

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) ;
     if ( s->ctrl.traveling_pseudo ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  TravelingPseudo(Sample *sample) : s(sample) { s->ctrl.traveling_pseudo = false; }
};
#endif

// Local Variables:
// mode: c++
// End:
