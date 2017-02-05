/*                            -*- C++ -*-
 * Copyright (C) 2012 Felix Salfelder
 * Author: Felix Salfelder <felix@salfelder.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * sensitivity analysis
 */
#include "u_sim_data.h"
#include "u_status.h"
#include "u_parameter.h"
#include "u_prblst.h"
#include "s__.h"
#include "io_matrix.h"
#include "e_node.h"
#include "e_aux.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class SENS : public SIM {
public:
  void	do_it(CS&, CARD_LIST*);

  explicit SENS():
    SIM(),
    _start(),
    _stop(),
    _step_in(),
    _step(0.),
    _linswp(false),
    _prevopppoint(false),
    _stepmode(ONE_PT),
    _dump_matrix(0)
  {}

  ~SENS() {}
private:
  explicit SENS(const SENS&):SIM() {unreachable(); incomplete();}
  void	sweep();
  void	first();
  bool	next();
  bool	next_output();
  bool	next_freq();
  void	solve();
  void	clear();
  void	setup(CS&);
private:
  PARAMETER<double> _start;	// sweep start frequency
  PARAMETER<double> _stop;	// sweep stop frequency
  PARAMETER<double> _step_in;	// step size, as input
  double _step;			// printed step size
  bool	_linswp;		// flag: use linear sweep (vs log sweep)
  bool	_prevopppoint;  	// flag: use previous op point
  enum {ONE_PT, LIN_STEP, LIN_PTS, TIMES, OCTAVE, DECADE} _stepmode;
  bool _dump_matrix; // dump matrix after ac
  struct output_t{
    string label;
    CKT_BASE* brh[2];
  };
  std::vector<output_t> _output;
  std::vector<output_t>::iterator _output_iter;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void SENS::do_it(CS& Cmd, CARD_LIST* Scope)
{
  _scope = Scope;
  _sim->set_command_sens();
  reset_timers();
  ::status.ac.reset().start();

  _sim->init();
  _sim->alloc_vectors();
  _sim->_acx.reallocate();
  _sim->_acx.set_min_pivot(OPT::pivtol);

  setup(Cmd);
  ::status.set_up.stop();
  switch (ENV::run_mode) {
  case rPRE_MAIN:	unreachable();	break;
  case rPIPE:		untested();
  case rBATCH:
  case rINTERACTIVE:
  case rSCRIPT:		sweep();	break;
  case rPRESET:		/*nothing*/	break;
  }
  _sim->_acx.unallocate();
  _sim->unalloc_vectors();

  // _output.clear();

  ::status.ac.stop();
  ::status.total.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
static int needslinfix;	// flag: lin option needs patch later (spice compat)
/*--------------------------------------------------------------------------*/
void SENS::setup(CS& Cmd)
{
  _out = IO::mstdout;
  _out.reset(); //BUG// don't know why this is needed

  bool ploton = IO::plotset  &&  plotlist().size() > 0;

// FIXME: freq-range?
//
  std::string output;
  bool newprobes = true;

  unsigned here = Cmd.cursor();
  do{
    // FIXME: can we use the probe parser?!
    if( Cmd >> "v" ){
      output_t t;
      t.brh[1] = 0;
      trace1("SENS::setup, have v", Cmd.tail());

      int paren = Cmd.skip1b('(');
      Cmd >> output;
      t.label = "v(" + output;

      CKT_NODE* node = dynamic_cast<CKT_NODE*>((*_scope).node(output));
      if(node)
        t.brh[0] = node;
      else{
        continue;
      }

      trace2("SENS::setup, have", output, Cmd.tail());
      if(!Cmd.match1(')')){
        output = Cmd.ctos(")");
        trace1("SENS:skip1b:setup, have2", output);
        t.label += "," + output;
        node = dynamic_cast<CKT_NODE*>((*_scope).node(output));
        if(node)
          t.brh[1] = node;
        else{
          Cmd.warn(bWARNING, "probelist: what's this?");
        }

      }else{
        trace0("SENS::setup no comma");
      }

      paren -= Cmd.skip1b(')');
      if (paren != 0) {untested();
        Cmd.warn(bWARNING, "need )");
      }else if (output.empty()) {untested();
        Cmd.warn(bWARNING, "probelist: what's this?");
      }else{
      }

      t.label+=")";
      if (newprobes){
        _output.clear();
        newprobes = false;
      }
      _output.push_back(t);
    }else if (Cmd.match1("'\"({") || Cmd.is_float()) { untested();
      Cmd >> _start;
      trace1("got start", _start);
      if (Cmd.match1("'\"({") || Cmd.is_float()) { untested();
        Cmd >> _stop;
      }else{ untested();
        _stop = _start;
      }
      if (Cmd.match1("'\"({") || Cmd.is_float()) { untested();
        _stepmode = LIN_STEP;
        Cmd >> _step_in;
      }else{
      }
    }else{
      incomplete();
    }
    //try{
    //  _output.add_list(Cmd);
    //}
    //catch(Exception_Cant_Find)
    //{}
    ONE_OF
      || (Get(Cmd, "*",		  &_step_in) && (_stepmode = TIMES))
      || (Get(Cmd, "+",		  &_step_in) && (_stepmode = LIN_STEP))
      || (Get(Cmd, "by",	  &_step_in) && (_stepmode = LIN_STEP))
      || (Get(Cmd, "step",	  &_step_in) && (_stepmode = LIN_STEP))
      || (Get(Cmd, "d{ecade}",	  &_step_in) && (_stepmode = DECADE))
      || (Get(Cmd, "ti{mes}",	  &_step_in) && (_stepmode = TIMES))
      || (Get(Cmd, "lin",	  &_step_in) && (_stepmode = LIN_PTS))
      || (Get(Cmd, "o{ctave}",	  &_step_in) && (_stepmode = OCTAVE))
      || Get(Cmd, "sta{rt}",	  &_start)
      || Get(Cmd, "sto{p}",	  &_stop)
      || Get(Cmd, "dm",	          &_dump_matrix)
      || _out.outset(Cmd);
      ;
  }while (Cmd.more() && !Cmd.stuck(&here));
  Cmd.check(bWARNING, "what's this??");

  _start.e_val(0., _scope);
  trace1("eval start", _start);
  _stop.e_val(0., _scope);
  _step_in.e_val(0., _scope);
  _step = _step_in;

  if (needslinfix) {untested();		// LIN option is # of points.
    assert(_step >= 2);			// Must compute step after 
    _step=(_stop-_start)/(_step-1.);	// reading start and stop,
    needslinfix = false;		// but step must be read first
  }else{			// for Spice compatibility
  }		
  if (_step==0.) { untested();
    _step = _stop - _start;
    _linswp = true;
  }else{
  }

  IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
  initio(_out);

  if(!_output.size()){
    throw(Exception("no output specified"));
  }
  trace2("eval done", _start, _stop);
}
/*--------------------------------------------------------------------------*/
void SENS::solve()
{
  _sim->_acx.zero();
  std::fill_n(_sim->_ac, _sim->_total_nodes+1, 0.);

  ::status.load.start();
  _sim->count_iterations(iTOTAL);
  CARD_LIST::card_list.do_ac();
  CARD_LIST::card_list.ac_load();
  ::status.load.stop();

  if (_dump_matrix){
    _out.setfloatwidth(0,0);
    _out << _sim->_acx << "\n" ;
  }
  ::status.lud.start();
  _sim->_acx.lu_decomp();
  ::status.lud.stop();

  // is this necessary?
  ::status.back.start();
  _sim->_acx.fbsub(_sim->_ac);
  ::status.back.stop();

  CKT_NODE* np = dynamic_cast<CKT_NODE*>((*_output_iter).brh[0]);
  assert(np);
  CKT_NODE* nn = dynamic_cast<CKT_NODE*>((*_output_iter).brh[1]);
  if(!nn) nn = &ground_node;

  trace2("set_sens_port", node_t(np).m_(), node_t(nn).m_());
  for(unsigned i=0; i<_sim->_total_nodes+1;i++){
    _sim->_sens[i] = 0.;
  }

  set_sens_port(node_t(np), node_t(nn));

  if(_dump_matrix){
    _out.setfloatwidth(0,0);
    _out << "sens_in\n";
    for(unsigned i=0; i<_sim->_total_nodes+1;i++){
      _out << _sim->_sens[i];
    }
    _out << "\n";
  }
  ::status.back.start();
  CKT_BASE::_sim->_acx.fbsubt(CKT_BASE::_sim->_sens);
  ::status.back.stop();

  if(_dump_matrix){
    _out.setfloatwidth(0,0);
    _out << "sens_out\n";
    for(unsigned i=0; i<_sim->_total_nodes+1;i++){
      _out << _sim->_sens[i];
    }
    _out << "\n";
  }

  CARD_LIST::card_list.do_sens();
}
/*--------------------------------------------------------------------------*/
void SENS::sweep()
{
  int width = std::min(OPT::numdgt+5, BIGBUFLEN-10);
  char format[20];
  //sprintf(format, "%%c%%-%u.%us", width, width);
  sprintf(format, "%%c%%-%us", width);

  _out.form(format, '*', "param");
  sprintf(format, "%%-%us", width);
  head(_start, _stop, "@freq");
  _sim->_jomega = 0; // COMPLEX(0., _sim->_freq * M_TWO_PI);
//  _sim->_jomega = COMPLEX(0., M_TWO_PI);
//
  trace2("eval done", _start, _stop);
  first();
  CARD_LIST::card_list.ac_begin();

  do { untested();
    do { untested();
      _sim->_jomega = COMPLEX(0., _sim->_freq * M_TWO_PI);
      solve();
      _out.form(format, (*_output_iter).label.c_str() );
      outdata(_sim->_freq, ofPRINT);
    } while (next_output());
  } while (next_freq());
}
/*--------------------------------------------------------------------------*/
void SENS::first()
{
  _sim->_freq = _start;
  assert(_output.size());
  _output_iter = _output.begin();
}
/*--------------------------------------------------------------------------*/
bool SENS::next_output()
{
  _output_iter++;
  return _output_iter != _output.end();
}
/*--------------------------------------------------------------------------*/
bool SENS::next_freq()
{ untested();
  _output_iter = _output.begin();
  double realstop = (_linswp)
    ? _stop - _step/100.
    : _stop / pow(_step,.01);
  if (!in_order(double(_start), _sim->_freq, realstop)) { untested();
    return false;
  }else{ untested();
  }

  _sim->_freq = (_linswp)
    ? _sim->_freq + _step
    : _sim->_freq * _step;
  if (in_order(_sim->_freq, double(_start), double(_stop))) { untested();
    return false;
  }else{
    return true;
  }
}
/*--------------------------------------------------------------------------*/
static SENS p1;
static DISPATCHER<CMD>::INSTALL d1(&command_dispatcher, "sens", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:et:
