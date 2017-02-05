/*                             -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 * Author: Albert Davis <aldavis@gnu.org>
 *         Felix Salfelder <felix@salfelder.org>
 *
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
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
 * alternative dc analysis
 */
#include "globals.h"
#include "u_status.h"
#include "u_prblst.h"
#include "u_cardst.h"
#include "e_elemnt.h"
#include "s__.h"
/*--------------------------------------------------------------------------*/
// bdaaac90de9b6a527b8d1c502e17fdf3bb2ab9bc
namespace {
/*--------------------------------------------------------------------------*/
class DCOP : public SIM {
protected:
  enum STEP_CAUSE {
    scUSER      =  1,	/* user requested				*/
    scSKIP      =  3,	/* effect of "skip" parameter			*/
    scITER_R    =  4,	/* iter count exceeds itl4 (reducing)		*/
    scITER_A    =  5,	/* iter count exceeds itl3 (holding)		*/
    scTE        =  6,	/* truncation error, or device stuff		*/
    scINITIAL   =  9,	/* initial guess				*/
    scREJECT    = 10,	/* rejected previous time step			*/
    scZERO      = 11,	/* fixed zero time step				*/
    scSMALL     = 12,	/* time step too small				*/
    scNO_ADVANCE= 13,	/* after all that it still didn't advance	*/
    scLAST      = 15 	/* last step */
  };
  void	set_step_cause(STEP_CAUSE C) {::status.control = C;}
  STEP_CAUSE step_cause()const {return STEP_CAUSE(::status.control);}
public:
  void	finish();
protected:
  void	fix_args(int);
  void	options(CS&, int);
  std::string label()const{untested(); return "dc";}
private:
  void	sweep();
  void	sweep_recursive(int);
  void	first(int);
  bool	next(int);
  explicit DCOP(const DCOP&): SIM() {unreachable(); incomplete();}
protected:
  explicit DCOP();
  ~DCOP() {}
  
protected:
  enum {DCNEST = 4};
  int _n_sweeps;
  PARAMETER<double> _start[DCNEST];
  PARAMETER<double> _stop[DCNEST];
  PARAMETER<double> _step_in[DCNEST];
  double _step[DCNEST];
  double _val_by_user_request[DCNEST];
  double _sweepdamp[DCNEST];
  bool _linswp[DCNEST];
  double* (_sweepval[DCNEST]);	/* pointer to thing to sweep, dc command */
  ELEMENT* (_zap[DCNEST]);	/* to branch to zap, for re-expand */
  CARDSTASH _stash[DCNEST];	/* store std values of elements being swept */
  bool _loop[DCNEST];		/* flag: do it again backwards */
  bool _reverse_in[DCNEST];	/* flag: sweep backwards, input */
  bool _reverse[DCNEST];	/* flag: sweep backwards, working */
  bool _cont;			/* flag: continue from previous run */
  TRACE _trace;			/* enum: show extended diagnostics */
  bool _converged;
  bool _ever_converged;         /* don't try to step back otherwise... */
  enum {ONE_PT, LIN_STEP, LIN_PTS, TIMES, OCTAVE, DECADE} _stepmode[DCNEST];
};
/*--------------------------------------------------------------------------*/
class DC : public DCOP {
public:
  explicit DC(): DCOP() {}
  ~DC() {}
  void	do_it(CS&, CARD_LIST*);
private:
  void	setup(CS&);
  explicit DC(const DC&): DCOP() {unreachable(); incomplete();}
};
/*--------------------------------------------------------------------------*/
class OP : public DCOP {
public:
  explicit OP(): DCOP() {}
  ~OP() {}
  void	do_it(CS&, CARD_LIST*);
private:
  void	setup(CS&);
  explicit OP(const OP&): DCOP() {unreachable(); incomplete();}
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DC::do_it(CS& Cmd, CARD_LIST* Scope)
{ untested();
  _scope = Scope;
  _sim->_time0 = 0.;
  _sim->set_command_dc();
  _sim->_phase = p_INIT_DC;
  ::status.dc.reset().start();
  command_base(Cmd);
  _sim->_has_op = s_DC;
  _scope = NULL;
  ::status.dc.stop();
}
/*--------------------------------------------------------------------------*/
void OP::do_it(CS& Cmd, CARD_LIST* Scope)
{ untested();
  _scope = Scope;
  _sim->_time0 = 0.;
  _sim->set_command_op();
  _sim->_phase = p_INIT_DC;
  ::status.op.reset().start();
  command_base(Cmd);
  _sim->_has_op = s_OP;
  _scope = NULL;
  ::status.op.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DCOP::DCOP()
  :SIM(),
   _n_sweeps(1),
   _cont(false),
   _trace(tNONE),
   _converged(0)
{ untested();
  for (int ii = 0; ii < DCNEST; ++ii) { untested();
    _loop[ii] = false;
    _reverse_in[ii] = false;
    _reverse[ii] = false;
    _step[ii]=0.;
    _linswp[ii]=true;
    _sweepval[ii]=&_sim->_genout;
    _zap[ii]=NULL; 
    _stepmode[ii] = ONE_PT;
  }
  
  //BUG// in SIM.  should be initialized there.
  //_sim->_genout=0.;
  _out=IO::mstdout;
}
/*--------------------------------------------------------------------------*/
void DCOP::finish(void)
{ untested();
  for (int ii = 0;  ii < _n_sweeps;  ++ii) { untested();
    if (_zap[ii]) { // component
      _stash[ii].restore();
      _zap[ii]->dec_probes();
      _zap[ii]->precalc_first();
      _zap[ii]->precalc_last();
      _zap[ii] = NULL;
    }else{ untested();
    }
  }
}
/*--------------------------------------------------------------------------*/
void OP::setup(CS& Cmd)
{ untested();
  _sim->_temp_c = OPT::temp_c;
  _cont = false;
  _trace = tNONE;
  _out = IO::mstdout;
  _out.reset(); //BUG// don't know why this is needed */
  bool ploton = IO::plotset  &&  plotlist().size() > 0;

  _zap[0] = NULL;
  _sweepval[0] = &(_sim->_temp_c);

  if (Cmd.match1("'\"({") || Cmd.is_float()) { untested();
    Cmd >> _start[0];
    if (Cmd.match1("'\"({") || Cmd.is_float()) { untested();
      Cmd >> _stop[0];
    }else{ untested();
      _stop[0] = _start[0];
    }
  }else{ untested();
  }
  
  _step[0] = 0.;
  _sim->_genout = 0.;

  options(Cmd,0);
//  _sim->_temp_c = temp_c_in;

  _n_sweeps = 1;
  Cmd.check(bWARNING, "what's this?");
  _sim->_freq = 0;

  IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
  initio(_out);

  _start[0].e_val(OPT::temp_c, _scope);
  fix_args(0);
}
/*--------------------------------------------------------------------------*/
void DC::setup(CS& Cmd)
{ untested();
  _sim->_temp_c = OPT::temp_c;
  _cont = false;
  _trace = tNONE;
  _out = IO::mstdout;
  _out.reset(); //BUG// don't know why this is needed */
  bool ploton = IO::plotset  &&  plotlist().size() > 0;

  if (Cmd.more()) { untested();
    for (_n_sweeps = 0; Cmd.more() && _n_sweeps < DCNEST; ++_n_sweeps) { untested();
      CARD_LIST::fat_iterator ci = findbranch(Cmd, &CARD_LIST::card_list);
      if (!ci.is_end()) {			// sweep a component
	if (ELEMENT* c = dynamic_cast<ELEMENT*>(*ci)) { untested();
	  _zap[_n_sweeps] = c;
	}else{ untested();
	  throw Exception("dc/op: can't sweep " + (**ci).long_label() + '\n');
	}
      }else if (Cmd.is_float()) {		// sweep the generator
	_zap[_n_sweeps] = NULL;
      }else{ untested();
	// leave as it was .. repeat Cmd with no args
      }
      
      if (Cmd.match1("'\"({") || Cmd.is_float()) {	// set up parameters
	_start[_n_sweeps] = "NA";
	_stop[_n_sweeps] = "NA";
	Cmd >> _start[_n_sweeps] >> _stop[_n_sweeps];
	_step[_n_sweeps] = 0.;
      }else{ untested();
	// leave it as it was .. repeat Cmd with no args
      }
      
      _sim->_genout = 0.;
      options(Cmd,_n_sweeps);
//      _sim->_temp_c = temp_c_in;
    }
  }else{ 
  }
  Cmd.check(bWARNING, "what's this?");

  IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
  initio(_out);

  assert(_n_sweeps > 0);
  for (int ii = 0;  ii < _n_sweeps;  ++ii) { untested();
    _start[ii].e_val(0., _scope);
    fix_args(ii);

    if (_zap[ii]) { // component
      _stash[ii] = _zap[ii];			// stash the std value
      _zap[ii]->inc_probes();			// we need to keep track of it
      _zap[ii]->set_value(_zap[ii]->value(),0);	// zap out extensions
      _zap[ii]->set_constant(false);		// so it will be updated
      _sweepval[ii] = _zap[ii]->set__value();	// point to value to patch
    }else{ // generator
      _sweepval[ii] = &_sim->_genout;			// point to value to patch
    }
  }
  _sim->_freq = 0;
}
/*--------------------------------------------------------------------------*/
void DCOP::fix_args(int Nest)
{ untested();
  _stop[Nest].e_val(_start[Nest], _scope);
  _step_in[Nest].e_val(0., _scope);
  _step[Nest] = _step_in[Nest];
  
  switch (_stepmode[Nest]) { untested();
  case ONE_PT:
  case LIN_STEP:
    _linswp[Nest] = true;
    break;
  case LIN_PTS:untested();
    if (_step[Nest] <= 2.) { untested();
      _step[Nest] = 2.;
    }else{ untested();
    }
    _linswp[Nest] = true;
    break;
  case TIMES:
    if (_step[Nest] == 0.  &&  _start[Nest] != 0.) { untested();
      _step[Nest] = _stop[Nest] / _start[Nest];
    }else{ untested();
    }
    _linswp[Nest] = false;
    break;
  case OCTAVE:untested();
    if (_step[Nest] == 0.) {untested();
      _step[Nest] = 1.;
    }else{untested();
    }
    _step[Nest] = pow(2.00000001, 1./_step[Nest]);
    _linswp[Nest] = false;
    break;
  case DECADE:
    if (_step[Nest] == 0.) {untested();
      _step[Nest] = 1.;
    }else{ untested();
    }
    _step[Nest] = pow(10., 1./_step[Nest]);
    _linswp[Nest] = false;
    break;
  };
  
  if (_step[Nest] == 0.) {	// prohibit log sweep from 0
    _step[Nest] = _stop[Nest] - _start[Nest];
    _linswp[Nest] = true;
  }else{ untested();
  }
}
/*--------------------------------------------------------------------------*/
void DCOP::options(CS& Cmd, int Nest)
{ untested();
  _sim->_uic = _loop[Nest] = _reverse_in[Nest] = false;
  unsigned here = Cmd.cursor();
  do{ untested();
    ONE_OF
      || (Cmd.match1("'\"({")	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
      || (Cmd.is_float()	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
      || (Get(Cmd, "*",		  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
      || (Get(Cmd, "+",		  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
      || (Get(Cmd, "by",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
      || (Get(Cmd, "step",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
      || (Get(Cmd, "d{ecade}",	  &_step_in[Nest]) && (_stepmode[Nest] = DECADE))
      || (Get(Cmd, "ti{mes}",	  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
      || (Get(Cmd, "lin",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_PTS))
      || (Get(Cmd, "o{ctave}",	  &_step_in[Nest]) && (_stepmode[Nest] = OCTAVE))
      || Get(Cmd, "c{ontinue}",   &_cont)
      || Get(Cmd, "uic",	  &_sim->_uic)
      || Get(Cmd, "dt{emp}",	  &(_sim->_temp_c),   mOFFSET, OPT::temp_c)
      || Get(Cmd, "lo{op}", 	  &_loop[Nest])
      || Get(Cmd, "re{verse}",	  &_reverse_in[Nest])
      || Get(Cmd, "te{mperature}",&(_sim->_temp_c))
      || (Cmd.umatch("tr{ace} {=}") &&
	  (ONE_OF
	   || Set(Cmd, "n{one}",      &_trace, tNONE)
	   || Set(Cmd, "o{ff}",       &_trace, tNONE)
	   || Set(Cmd, "w{arnings}",  &_trace, tUNDER)
	   || Set(Cmd, "a{lltime}",   &_trace, tALLTIME)
	   || Set(Cmd, "r{ejected}",  &_trace, tREJECTED)
	   || Set(Cmd, "i{terations}",&_trace, tITERATION)
	   || Set(Cmd, "v{erbose}",   &_trace, tVERBOSE)
	   || Cmd.warn(bWARNING, 
		       "need none, off, warnings, iterations, verbose")
	   )
	  )
      || outset(Cmd,&_out);
  }while (Cmd.more() && !Cmd.stuck(&here));
}
/*--------------------------------------------------------------------------*/
void DCOP::sweep()
{ untested();
  head(_start[0], _stop[0], " ");
  _sim->_bypass_ok = false;
  _sim->set_inc_mode_bad();
  if (_cont) { untested();
    _sim->restore_voltages();
    CARD_LIST::card_list.tr_restore();
  }else{ untested();
    _sim->clear_limit();
    CARD_LIST::card_list.tr_begin();
  }

  set_step_cause(scUSER);
  _converged = false;
  _ever_converged = false;
  for (int ii = 0; ii < _n_sweeps; ++ii) { untested();
    if (!_zap[ii]) { untested();
    }else if (_zap[ii]->is_constant()) { untested();
      incomplete();
//      CARD_LIST::card_list.q_hack(_zap[ii]);
    }else{ untested();
    }
  }
  sweep_recursive(_n_sweeps-1);
  _sim->pop_voltages();
  _sim->_has_op = _sim->_mode;
  _sim->keep_voltages();
}
/*--------------------------------------------------------------------------*/
static double mul(double a,double b){ return a*b; }
static double sub(double a,double b){ return a-b; }
static double div(double a,double b){ return a/b; }
static double add(double a,double b){ return a+b; }
static bool ge(double a,double b){ return a>=b; }
static bool le(double a,double b){ return a<=b; }
/*--------------------------------------------------------------------------*/
void DCOP::sweep_recursive(int Nest)
{ untested();
  static unsigned extra_steps;
  assert(Nest >= 0);
  assert(Nest < DCNEST);

  OPT::ITL itl = OPT::DCBIAS;

  first(Nest);

  double (*step)(double a, double b) = add;
  double (*back)(double a, double b) = sub;
  if (!_linswp[Nest]) { untested();
    step = mul;
    back = div;
  }
  if (_reverse[Nest]) { untested();
    std::swap(step,back);
  }
  
  trace3("DCOP::sweep_recursive", Nest, *(_sweepval[Nest]), _step[Nest]);

  bool firstloop=true;
  do { untested();
//     _sim->_temp_c = temp_c_in;
    if (Nest) { untested();
      sweep_recursive(Nest-1);
      if(_converged || !_ever_converged){ untested();
	//if(!firstloop) 
	_sim->pop_voltages();
	_sim->restore_voltages();
        if(firstloop) { untested();
	  _sim->_has_op = _sim->_mode;
	  _sim->push_voltages();
	  _sim->keep_voltages();
	}
      } else { untested();
	trace2("not converged II", Nest, *(_sweepval[Nest]));
	// step back...
	incomplete();
      }
    }else{ // leaf
      if (_sim->command_is_op()) { untested();
	CARD_LIST::card_list.precalc_last();
      }else{ untested();
      }
      for (int ii = 0;  ii < _n_sweeps;  ++ii) { untested();
	if (_zap[ii]) { untested();
	//  _zap[ii]->do_tr(); // uf?
	}
      }
      _converged = solve_with_homotopy(itl,_trace);
      _ever_converged |= _converged;
      ++::status.hidden_steps;
      bool printnow =
        (_trace >= tREJECTED)
//	|| firstloop //uf?
        || (_converged && ((_trace >= tALLTIME)
        || (step_cause() == scUSER )));
      trace4("outdata?", _trace, _converged, printnow, _ever_converged);
      if (!printnow) { untested();
	++extra_steps;
	if(extra_steps > 100){ untested();
	  throw Exception("dc stepping did not succeed");
	}
      }
      if (!_converged) { untested();
	trace3("DCOP::sweep_recursive noconv", Nest, *_sweepval[Nest], _ever_converged);
	error(bWARNING, "did not converge\n");
        _sim->restore_voltages();
      }else if (!_ever_converged) { incomplete();
        if (firstloop) { untested();
	  _sim->push_voltages();
	  _sim->keep_voltages();
	}
      }else{ untested();
        ::status.accept.start();
        _sim->set_limit();
        CARD_LIST::card_list.tr_accept();
        ::status.accept.stop();
        _sim->keep_voltages();
        if (firstloop) { untested();
	  _sim->_has_op = _sim->_mode;
	  _sim->push_voltages();
	  _sim->keep_voltages();
	}
      }
      if (printnow) { untested();
	extra_steps = 0;
	fixzero(_sweepval[Nest], _step[Nest]); // hack
	if (_converged){ untested();
	  outdata(*_sweepval[Nest], ofPRINT | ofSTORE);
	} else { untested();
	  outdata(- *_sweepval[Nest], ofPRINT | ofSTORE);
	}
	::status.hidden_steps = 0;
      }else{ untested();
      }

      if (!_converged && firstloop && Nest) { untested();
	trace1("didnt converge in first", Nest);
	return;
      }
#if 0
      ::status.accept.start();
      _sim->set_limit();
      CARD_LIST::card_list.tr_accept();
      ::status.accept.stop();
      _sim->_has_op = _sim->_mode;
      outdata(*_sweepval[Nest], ofPRINT | ofSTORE | ofKEEP);
#endif
      itl = OPT::DCXFER;

    }

    if(firstloop){ untested();
      if(step_cause() != scUSER){ untested();
	trace1("firststep nouser", Nest);
	return;
      }
    } else { untested();
      // UGLY. next may have changed _reverse[Nest]
      step = add;
      back = sub;
      if (!_linswp[Nest]) { untested();
	step=mul;
	back=div;
      }
      if (_reverse[Nest]) { untested();
	std::swap(step,back);
      }
      // /UGLY

    }
    if ((firstloop || _converged) && step_cause() == scUSER) { untested();
      _val_by_user_request[Nest] = step(_val_by_user_request[Nest], _step[Nest]);
      trace2("ordered next step loop", Nest, _val_by_user_request[Nest]);
    }
    firstloop = false;
  } while (next(Nest));

  // _sim->pop_voltages();
}
/*--------------------------------------------------------------------------*/
void DCOP::first(int Nest)
{ untested();
  assert(Nest >= 0);
  assert(Nest < DCNEST);
  assert(_start);
  assert(_sweepval);
  assert(_sweepval[Nest]);

  if (ELEMENT* c = dynamic_cast<ELEMENT*>(_zap[Nest])) { untested();
    c->set_constant(false);
    // because of extra precalc_last could set constant to true
    // will be obsolete, once pointer hack is fixed
  }else{ untested();
    // not needed if not sweeping an element
  }

  *_sweepval[Nest] = _start[Nest];
  if(_converged){ untested();
    trace1("BUG", Nest);
//    set_step_cause(scUSER);
  }else{ untested();
  }
  assert(step_cause());

  _val_by_user_request[Nest] = _start[Nest];
  _sweepdamp[Nest] = 1;
  _reverse[Nest] = false;
  if (_reverse_in[Nest]) { untested();
    _converged = true; //?
    double (*step)(double a, double b) = add;
    double (*back)(double a, double b) = sub;
    if (!_linswp[Nest]) { untested();
      step=mul;
      back=div;
    }
    while (next(Nest)) { untested();
      _val_by_user_request[Nest] = step(*(_sweepval[Nest]), _step[Nest]);
    }
    _val_by_user_request[Nest] = back(*(_sweepval[Nest]), _step[Nest]);
    _reverse[Nest] = true;
    next(Nest);
    _converged = false;
  }else{ untested();
  }
  _sim->_phase = p_INIT_DC;
}
/*--------------------------------------------------------------------------*/
bool DCOP::next(int Nest)
{ itested();
  double sweepval = NOT_VALID;
  trace2("next", Nest, *(_sweepval[Nest]));
  bool ok = false;
  double nothing = 0;
  double (*step)(double a, double b) = add;
  double (*back)(double a, double b) = sub;
  bool (*further)(double a, double b) = ge;
  double (*scale)(double a, double b) = mul;
  if (!_linswp[Nest]) { untested();
    step = mul;
    back = div;
    scale = pow;
    nothing = 1;
    scale = pow;
  }
  double fudge = scale(_step[Nest], 1e-6);
  if (_reverse[Nest]) { untested();
    trace2("next, reverse", *_sweepval[Nest], _val_by_user_request[Nest]);
    std::swap(step,back);
    if (_step[Nest] < 0) { untested();
      further = ge;
    }else{ itested();
      further = le;
    }
    fudge = scale(_step[Nest], -1e-6);
  } else if (_step[Nest] < 0) { untested();
    further = le;
  }
  if (_step[Nest] == nothing) { untested();
    ok = false;
    set_step_cause(scZERO);
  }else if (!_converged && _ever_converged) { untested();
    if (_sweepdamp[Nest]<OPT::dtmin) { untested();
      throw Exception("step too small (does not converge)");
    }else{ untested();
    }
    _sweepdamp[Nest] /= 2.;
    trace2("reducing step by", _sweepdamp[Nest], Nest);
    sweepval = back(*_sweepval[Nest], scale(_step[Nest],_sweepdamp[Nest]));
    trace2("next at", sweepval, _val_by_user_request[Nest]);
    ok = true;
    set_step_cause(scREJECT);
  }else{ untested();
    if (_sweepdamp[Nest] != 1) { untested();
      trace3("recovered from", _sweepdamp[Nest], Nest, sweepval);
      set_step_cause(scTE);
    }else{ untested();
    }
    _sweepdamp[Nest] *= 1.4;
    _sweepdamp[Nest] = std::min(_sweepdamp[Nest],1.);
    sweepval = step(*_sweepval[Nest], scale(_step[Nest],_sweepdamp[Nest]));
    fixzero(&sweepval, _step[Nest]);
    ok = in_order(back(_start[Nest],fudge), sweepval, step(_stop[Nest],fudge));
    if(ok){ untested();
    }else if (!_reverse[Nest] && !ok && _loop[Nest]) { untested();
      _reverse[Nest] = true;
      if (_step[Nest] < 0) { untested();
	further = le;
      }else{ untested();
	further = ge;
      }
      fudge = scale(_step[Nest], -.1);
      sweepval = back(*_sweepval[Nest], _step[Nest]);
      std::swap(step,back);
      ok = in_order(back(_start[Nest],fudge), sweepval, step(_stop[Nest],fudge));
      assert(ok);
      trace2("BUG?", sweepval, *_sweepval[Nest]);
      _val_by_user_request[Nest] = sweepval; // BUG: here?
    }else if(_reverse_in[Nest]){ untested();
      // hmm maybe _reverse_in && !_reverse?
      trace2("reverse end?", sweepval, *_sweepval[Nest]);
      *_sweepval[Nest] = sweepval;
    }else{ untested();
    }
  }

  double v = _val_by_user_request[Nest];
  if(!ok){ untested();
  }else if (further(step(sweepval, scale(_step[Nest],1e-6) ), v)) { untested();
    trace5("userstep at", v, sweepval, ok, _reverse[Nest], *_sweepval[Nest]);
    set_step_cause(scUSER); // here?!
    sweepval = v;
  }else{ untested();
  }

  _sim->_phase = p_DC_SWEEP;
//  *(_sweepval[Nest]) = sweepval; // ouch.
  if (ok) { untested();
    assert(sweepval != NOT_VALID);
    *(_sweepval[Nest]) = sweepval;
    return true;
  }else{ untested();
    trace3("not ok at", v, sweepval, *(_sweepval[Nest]));
    //assert(sweepval == NOT_VALID);
    return false;
  }
}
/*--------------------------------------------------------------------------*/
static DC p2;
static OP p4;
static DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "dc", &p2);
static DISPATCHER<CMD>::INSTALL d4(&command_dispatcher, "op", &p4);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
