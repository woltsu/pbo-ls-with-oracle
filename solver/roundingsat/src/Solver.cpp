/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2020, Jo Devriendt
Copyright (c) 2020, Stephan Gocht

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

#include "Solver.hpp"
#include <cmath>
#include "Constr.hpp"
#include "aux.hpp"
#include "globals.hpp"
#include <algorithm>
#include <random>

namespace rs {

// ---------------------------------------------------------------------
// Initialization

Solver::Solver() : order_heap(activity) {
  ca.capacity(1024 * 1024);  // 4MB
}

void Solver::setNbVars(long long nvars, bool orig) {
  assert(nvars > 0);
  assert(nvars < INF);
  if (nvars <= n) return;
  aux::resizeIntMap(_adj, adj, nvars, resize_factor, {});
  aux::resizeIntMap(_Level, Level, nvars, resize_factor, INF);
  Pos.resize(nvars + 1, INF);
  Reason.resize(nvars + 1, CRef_Undef);
  activity.resize(nvars + 1, 1 / actLimitV);
  phase.resize(nvars + 1);
  cePools.resize(nvars + 1);
  order_heap.resize(nvars + 1);
  for (Var v = n + 1; v <= nvars; ++v) phase[v] = -v, order_heap.insert(v);
  // if (lpSolver) lpSolver->setNbVariables(nvars + 1); // Currently, LP solver only reasons on formula constraints
  n = nvars;
  if (orig) {
    orig_n = n;
    stats.NORIGVARS = n;
  }
}

void Solver::init() {
  if (!options.proofLog.get().empty()) logger = std::make_shared<Logger>(options.proofLog.get());
  cePools.initializeLogging(logger);
}

void Solver::initLP([[maybe_unused]] const CeArb objective) {
#if WITHSOPLEX
  if (options.lpPivotRatio.get() == 0) return;
  bool pureCNF = objective->vars.size() == 0;
  for (CRef cr : constraints) {
    if (!pureCNF) break;
    pureCNF = (ca[cr].degree() == 1);
  }
  if (pureCNF) return;
  lpSolver = std::make_shared<LpSolver>(*this, objective);
#else
  return;
#endif  // WITHSOPLEX
}

// ---------------------------------------------------------------------
// VSIDS

void Solver::vDecayActivity() { v_vsids_inc *= (1 / options.varDecay.get()); }
void Solver::vBumpActivity(Var v) {
  assert(v > 0);
  if ((activity[v] += v_vsids_inc) > actLimitV) {  // Rescale
    for (Var x = 1; x <= n; ++x) {
      activity[x] /= actLimitV;
      activity[x] /= actLimitV;
    }
    v_vsids_inc /= actLimitV;
    v_vsids_inc /= actLimitV;
  }
  // Update heap with respect to new activity:
  if (order_heap.inHeap(v)) order_heap.percolateUp(v);
}

void Solver::cDecayActivity() { c_vsids_inc *= (1 / options.clauseDecay.get()); }
void Solver::cBumpActivity(Constr& c) {
  c.act += c_vsids_inc;
  if (c.act > actLimitC) {  // Rescale:
    for (size_t i = 0; i < constraints.size(); i++) {
      ca[constraints[i]].act /= actLimitC;
      ca[constraints[i]].act /= actLimitC;
    }
    c_vsids_inc /= actLimitC;
    c_vsids_inc /= actLimitC;
  }
}

// ---------------------------------------------------------------------
// Assignment manipulation

void Solver::uncheckedEnqueue(Lit p, CRef from) {
  assert(!isTrue(Level, p));
  assert(!isFalse(Level, p));
  assert(isUnknown(Pos, p));
  Var v = toVar(p);
  Reason[v] = from;
  if (decisionLevel() == 0) {
    Reason[v] = CRef_Undef;  // no need to keep track of reasons for unit literals
    if (logger) {
      Constr& C = ca[from];
      C.toExpanded(cePools)->logUnit(Level, Pos, v, stats);
      assert(logger->unitIDs.size() == trail.size() + 1);
    }
  }
  Level[p] = decisionLevel();
  Pos[v] = (int)trail.size();
  trail.push_back(p);
}

void Solver::undoOne() {
  assert(!trail.empty());
  ++stats.NTRAILPOPS;
  Lit l = trail.back();
  if (qhead == (int)trail.size()) {
    for (const Watch& w : adj[-l])
      if (w.idx >= INF) ca[w.cref].undoFalsified(w.idx);
    --qhead;
  }
  Var v = toVar(l);
  trail.pop_back();
  Level[l] = INF;
  Pos[v] = INF;
  phase[v] = l;
  if (!trail_lim.empty() && trail_lim.back() == (int)trail.size()) trail_lim.pop_back();
  order_heap.insert(v);
}

void Solver::backjumpTo(int level) {
  assert(level >= 0);
  while (decisionLevel() > level) undoOne();
}

void Solver::decide(Lit l) {
  //++stats.NDECIDE;
  trail_lim.push_back(trail.size());
  uncheckedEnqueue(l, CRef_Undef);
}

void Solver::propagate(Lit l, CRef reason) {
  assert(reason != CRef_Undef);
  ++stats.NPROP;
  uncheckedEnqueue(l, reason);
}

/**
 * Unit propagation with watched literals.
 * @post: all watches up to trail[qhead] have been propagated
 */
CeSuper Solver::runPropagation(bool onlyUnitPropagation) {
  CeSuper confl = processLearnedStack();
  if (confl) {
    return confl;
  }
  while (qhead < (int)trail.size()) {
    Lit p = trail[qhead++];
    std::vector<Watch>& ws = adj[-p];
    for (int it_ws = 0; it_ws < (int)ws.size(); ++it_ws) {
      int idx = ws[it_ws].idx;
      if (idx < 0 && isTrue(Level, idx + INF)) {
        assert(dynamic_cast<Clause*>(&(ca[ws[it_ws].cref])) != nullptr);
        continue;
      }  // blocked literal check
      CRef cr = ws[it_ws].cref;
      WatchStatus wstat = checkForPropagation(cr, ws[it_ws].idx, -p);
      if (wstat == WatchStatus::DROPWATCH)
        aux::swapErase(ws, it_ws--);
      else if (wstat == WatchStatus::CONFLICTING) {  // clean up current level and stop propagation
        ++stats.NTRAILPOPS;
        for (int i = 0; i <= it_ws; ++i) {
          const Watch& w = ws[i];
          if (w.idx >= INF) ca[w.cref].undoFalsified(w.idx);
        }
        --qhead;
        Constr& C = ca[cr];
        if (!C.isLocked()) {
          cBumpActivity(C);
          recomputeLBD(C);
        }

        stats.NENCFORMULA += C.getOrigin() == Origin::FORMULA;
        stats.NENCLEARNED += C.getOrigin() == Origin::LEARNED;
        stats.NENCBOUND += (C.getOrigin() == Origin::UPPERBOUND);
        stats.NLPENCGOMORY += C.getOrigin() == Origin::GOMORY;
        stats.NLPENCLEARNEDFARKAS += C.getOrigin() == Origin::LEARNEDFARKAS;
        stats.NLPENCFARKAS += C.getOrigin() == Origin::FARKAS;

        return C.toExpanded(cePools);
      }
    }
  }
  if (onlyUnitPropagation) return CeNull();
  if (lpSolver) {
    std::pair<LpStatus, CeSuper> lpResult =
        aux::timeCall<std::pair<LpStatus, CeSuper>>([&] { return lpSolver->checkFeasibility(); }, stats.LPTOTALTIME);
    assert((lpResult.first == LpStatus::INFEASIBLE) == (lpResult.second && lpResult.second->hasNegativeSlack(Level)));
    return lpResult.second;
  }
  return CeNull();
  ;
}

WatchStatus Solver::checkForPropagation(CRef cr, int& idx, Lit p) {
  assert(isFalse(Level, p));
  Constr& C = ca[cr];
  if (C.isMarkedForDelete()) return WatchStatus::DROPWATCH;
  ++stats.NWATCHLOOKUPS;
  return C.checkForPropagation(cr, idx, p, *this);
}

// ---------------------------------------------------------------------
// Conflict analysis

void Solver::recomputeLBD(Constr& C) {
  if (C.lbd() > 2) {  // constraints with LBD <= 2 won't have score recomputed
    assert(tmpSet.isEmpty());
    for (int i = 0; i < (int)C.size(); i++)
      if (isFalse(Level, C.lit(i))) tmpSet.add(Level[-C.lit(i)]);
    if (C.lbd() > tmpSet.size() + 1) C.setLBD(tmpSet.size());  // simulate Glucose
    tmpSet.clear();
  }
}

CeSuper getAnalysisCE(const CeSuper& conflict, int bitsOverflow, ConstrExpPools& cePools) {
  if (bitsOverflow == 0 || bitsOverflow > conflLimit128) {
    CeArb confl = cePools.takeArb();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit96) {
    Ce128 confl = cePools.take128();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit64) {
    Ce96 confl = cePools.take96();
    conflict->copyTo(confl);
    return confl;
  } else if (options.bitsOverflow.get() > conflLimit32) {
    Ce64 confl = cePools.take64();
    conflict->copyTo(confl);
    return confl;
  } else {
    Ce32 confl = cePools.take32();
    conflict->copyTo(confl);
    return confl;
  }
}

CeSuper Solver::analyze(CeSuper conflict) {
  if (logger) logger->logComment("Analyze", stats);
  assert(conflict->hasNegativeSlack(Level));
  stats.NADDEDLITERALS += conflict->vars.size();
  conflict->removeUnitsAndZeroes(Level, Pos);
  conflict->saturateAndFixOverflow(getLevel(), (bool)options.weakenFull, options.bitsOverflow.get(),
                                   options.bitsReduced.get(), 0);

  CeSuper confl = getAnalysisCE(conflict, options.bitsOverflow.get(), cePools);
  conflict->reset();

  assert(actSet.isEmpty());  // will hold the literals that need their activity bumped
  for (Var v : confl->vars) {
    if (options.bumpLits)
      actSet.add(confl->getLit(v));
    else if (!options.bumpOnlyFalse || isFalse(Level, confl->getLit(v)))
      actSet.add(v);
  }
  while (decisionLevel() > 0) {
    if (asynch_interrupt) throw asynchInterrupt;
    Lit l = trail.back();
    if (confl->hasLit(-l)) {
      assert(confl->hasNegativeSlack(Level));
      AssertionStatus status = confl->isAssertingBefore(Level, decisionLevel());
      if (status == AssertionStatus::ASSERTING)
        break;
      else if (status == AssertionStatus::FALSIFIED) {
        backjumpTo(decisionLevel() - 1);
        assert(confl->hasNegativeSlack(Level));
        continue;
      }
      assert(isPropagated(Reason, l));
      Constr& reasonC = ca[Reason[toVar(l)]];
      if (!reasonC.isLocked()) {
        cBumpActivity(reasonC);
        recomputeLBD(reasonC);
      }

      stats.NENCFORMULA += reasonC.getOrigin() == Origin::FORMULA;
      stats.NENCLEARNED += reasonC.getOrigin() == Origin::LEARNED;
      stats.NENCBOUND += (reasonC.getOrigin() == Origin::UPPERBOUND);
      stats.NLPENCGOMORY += reasonC.getOrigin() == Origin::GOMORY;
      stats.NLPENCLEARNEDFARKAS += reasonC.getOrigin() == Origin::LEARNEDFARKAS;
      stats.NLPENCFARKAS += reasonC.getOrigin() == Origin::FARKAS;
      ++stats.NRESOLVESTEPS;

      reasonC.resolveWith(confl, l, &actSet, *this);
    }
    undoOne();
  }
  assert(confl->hasNegativeSlack(Level));
  for (Lit l : actSet.keys)
    if (l != 0) vBumpActivity(toVar(l));
  actSet.clear();

  return confl;
}

// ---------------------------------------------------------------------
// Constraint management

CRef Solver::attachConstraint(CeSuper constraint, bool locked) {
  assert(constraint->isSortedInDecreasingCoefOrder());
  assert(constraint->isSaturated());
  assert(constraint->hasNoZeroes());
  assert(constraint->hasNoUnits(getLevel()));
  assert(!constraint->isTautology());
  assert(constraint->vars.size() > 0);
  assert(!constraint->hasNegativeSlack(getLevel()));
  assert(constraint->orig != Origin::UNKNOWN);

  CRef cr = constraint->toConstr(ca, locked, logger ? constraint->logProofLineWithInfo("Attach", stats) : ++crefID);
  Constr& C = ca[cr];
  C.initializeWatches(cr, *this);
  constraints.push_back(cr);

  bool learned = (C.getOrigin() == Origin::LEARNED || C.getOrigin() == Origin::LEARNEDFARKAS ||
                  C.getOrigin() == Origin::FARKAS || C.getOrigin() == Origin::GOMORY);
  if (learned) {
    stats.LEARNEDLENGTHSUM += C.size();
    stats.LEARNEDDEGREESUM += C.degree();
  } else {
    stats.EXTERNLENGTHSUM += C.size();
    stats.EXTERNDEGREESUM += C.degree();
  }
  if (C.degree() == 1) {
    stats.NCLAUSESLEARNED += learned;
    stats.NCLAUSESEXTERN += !learned;
  } else if (C.largestCoef() == 1) {
    stats.NCARDINALITIESLEARNED += learned;
    stats.NCARDINALITIESEXTERN += !learned;
  } else {
    stats.NGENERALSLEARNED += learned;
    stats.NGENERALSEXTERN += !learned;
  }

  stats.NCONSFORMULA += C.getOrigin() == Origin::FORMULA;
  stats.NCONSLEARNED += C.getOrigin() == Origin::LEARNED;
  stats.NCONSBOUND += (C.getOrigin() == Origin::UPPERBOUND);
  stats.NLPGOMORYCUTS += C.getOrigin() == Origin::GOMORY;
  stats.NLPLEARNEDFARKAS += C.getOrigin() == Origin::LEARNEDFARKAS;
  stats.NLPFARKAS += C.getOrigin() == Origin::FARKAS;

  return cr;
}

void Solver::learnConstraint(const CeSuper c, Origin orig) {
  assert(orig == Origin::LEARNED || orig == Origin::FARKAS || orig == Origin::LEARNEDFARKAS || orig == Origin::GOMORY);
  CeSuper learned = c->clone(cePools);
  learned->orig = orig;
  learned->saturateAndFixOverflow(getLevel(), (bool)options.weakenFull, options.bitsLearned.get(),
                                  options.bitsLearned.get(), 0);
  learnedStack.push_back(learned->toSimple());
}

// NOTE: backjumps to place where the constraint is not conflicting, as otherwise we might miss propagations
CeSuper Solver::processLearnedStack() {
  // loop back to front as the last constraint in the queue is a result of conflict analysis
  // and we want to first check this constraint to backjump.
  while (learnedStack.size() > 0) {
    CeSuper learned = learnedStack.back()->toExpanded(cePools);
    learnedStack.pop_back();
    learned->removeUnitsAndZeroes(Level, Pos);
    learned->sortInDecreasingCoefOrder();
    int assertionLevel = learned->getAssertionLevel(Level, Pos);
    if (assertionLevel < 0) {
      backjumpTo(0);
      return learned;
    }
    backjumpTo(assertionLevel);
    assert(!learned->hasNegativeSlack(Level));
    learned->heuristicWeakening(Level, Pos, stats);  // TODO: don't always weaken heuristically?
    learned->postProcess(Level, Pos, false, stats);
    assert(learned->isSaturated());
    if (learned->isTautology()) {
      continue;
    }
    CRef cr = attachConstraint(learned, false);
    Constr& C = ca[cr];
    if (assertionLevel < INF)
      recomputeLBD(C);
    else
      C.setLBD(C.size());  // the LBD of non-asserting constraints is undefined, so we take a safe upper bound
  }
  return CeNull();
}

std::pair<ID, ID> Solver::addInputConstraint(CeSuper ce) {
  assert(ce->orig == Origin::FORMULA || ce->orig == Origin::UPPERBOUND);
  assert(decisionLevel() == 0);
  ID input = ID_Undef;
  if (logger) input = ce->logAsInput();
  ce->postProcess(Level, Pos, true, stats);
  if (ce->isTautology()) {
    return {input, ID_Undef};  // already satisfied.
  }

  if (ce->hasNegativeSlack(Level)) {
    if (options.verbosity.get() > 0) puts("c Inconsistent input constraint");
    if (logger) ce->logInconsistency(Level, Pos, stats);
    assert(decisionLevel() == 0);
    return {input, ID_Unsat};
  }

  if (options.bitsInput.get() != 0 && !ce->largestCoefFitsIn(options.bitsInput.get())) {
    ce->toStreamAsOPB(std::cerr);
    quit::exit_ERROR({"Input constraint coefficient exceeds bit limit."});
  }

  CRef cr = attachConstraint(ce, true);
  CeSuper confl = aux::timeCall<CeSuper>([&] { return runPropagation(true); }, stats.PROPTIME);
  if (confl) {
    assert(confl->hasNegativeSlack(Level));
    if (options.verbosity.get() > 0) puts("c Input conflict");
    if (logger) confl->logInconsistency(Level, Pos, stats);
    assert(decisionLevel() == 0);
    return {input, ID_Unsat};
  }
  ID id = ca[cr].id;
  Origin orig = ca[cr].getOrigin();
  if (orig != Origin::FORMULA) {
    external[id] = cr;
  }
  if (lpSolver && (orig == Origin::FORMULA || orig == Origin::UPPERBOUND)) {
    lpSolver->addConstraint(cr, false, orig == Origin::UPPERBOUND, false);
  }
  return {input, id};
}

std::pair<ID, ID> Solver::addConstraint(const CeSuper c, Origin orig) {
  // NOTE: copy to temporary constraint guarantees original constraint is not changed and does not need logger
  CeSuper ce = c->clone(cePools);
  ce->orig = orig;
  std::pair<ID, ID> result = addInputConstraint(ce);
  return result;
}

std::pair<ID, ID> Solver::addConstraint(const ConstrSimpleSuper& c, Origin orig) {
  CeSuper ce = c.toExpanded(cePools);
  ce->orig = orig;
  std::pair<ID, ID> result = addInputConstraint(ce);
  return result;
}

std::pair<ID, ID> Solver::addUnitConstraint(Lit l, Origin orig) {
  return addConstraint(ConstrSimple32({{1, l}}, 1), orig);
}

void Solver::removeConstraint(Constr& C, [[maybe_unused]] bool override) {
  assert(override || !C.isLocked());
  assert(!C.isMarkedForDelete());
  assert(!external.count(C.id));
  C.markForDel();
  ca.wasted += C.getMemSize();
}

void Solver::dropExternal(ID id, bool erasable, bool forceDelete) {
  assert(erasable || !forceDelete);
  if (id == ID_Undef) return;
  auto old_it = external.find(id);
  assert(old_it != external.end());
  Constr& constr = ca[old_it->second];
  external.erase(old_it);
  constr.setLocked(!erasable);
  if (forceDelete) removeConstraint(constr);
}

// ---------------------------------------------------------------------
// Garbage collection

// We assume in the garbage collection method that reduceDB() is the
// only place where constraints are deleted.
void Solver::garbage_collect() {
  assert(decisionLevel() == 0);  // so we don't need to update the pointer of Reason<CRef>
  if (options.verbosity.get() > 0) puts("c GARBAGE COLLECT");

  ca.wasted = 0;
  ca.at = 0;
  std::unordered_map<uint32_t, CRef> crefmap;
  for (int i = 1; i < (int)constraints.size(); ++i) assert(constraints[i - 1].ofs < constraints[i].ofs);
  for (CRef& cr : constraints) {
    uint32_t offset = cr.ofs;
    size_t memSize = ca[cr].getMemSize();
    memmove(ca.memory + ca.at, ca.memory + cr.ofs, sizeof(uint32_t) * memSize);
    cr.ofs = ca.at;
    ca.at += memSize;
    crefmap[offset] = cr;
  }
#define update_ptr(cr) cr = crefmap[cr.ofs];
  for (Lit l = -n; l <= n; l++)
    for (size_t i = 0; i < adj[l].size(); i++) update_ptr(adj[l][i].cref);
  for (auto& ext : external) update_ptr(ext.second);
#undef update_ptr
}

// We assume in the garbage collection method that reduceDB() is the
// only place where constraints are removed from memory.
void Solver::reduceDB() {
  assert(decisionLevel() == 0);

  std::vector<CRef> learnts;
  learnts.reserve(constraints.size() / 2);

  size_t totalLearnts = 0;
  size_t promisingLearnts = 0;
  for (CRef& cr : constraints) {
    Constr& C = ca[cr];
    if (C.isMarkedForDelete() || external.count(C.id)) continue;
    BigVal eval = -C.degree();
    for (int j = 0; j < (int)C.size() && eval < 0; ++j)
      if (isUnit(Level, C.lit(j))) eval += C.coef(j);
    if (eval >= 0)
      removeConstraint(C, true);  // remove constraints satisfied at root
    else if (!C.isLocked()) {
      if (C.size() > 2 && C.lbd() > 2) learnts.push_back(cr);  // Keep all binary clauses and short LBDs
      if (C.size() <= 2 || C.lbd() <= 3) ++promisingLearnts;
      ++totalLearnts;
    }
  }

  if (promisingLearnts > totalLearnts / 2)
    nconfl_to_reduce += 10 * options.dbCleanInc.get();
  else
    nconfl_to_reduce += options.dbCleanInc.get();
  std::sort(learnts.begin(), learnts.end(), [&](CRef x, CRef y) {
    return ca[x].lbd() > ca[y].lbd() || (ca[x].lbd() == ca[y].lbd() && ca[x].act < ca[y].act);
  });
  for (size_t i = 0; i < std::min(totalLearnts / 2, learnts.size()); ++i) removeConstraint(ca[learnts[i]]);

  for (Lit l = -n; l <= n; ++l)
    for (int i = 0; i < (int)adj[l].size(); ++i) {
      if (ca[adj[l][i].cref].isMarkedForDelete()) aux::swapErase(adj[l], i--);
    }

  size_t i = 0;
  size_t j = 0;
  for (; i < constraints.size(); ++i) {
    Constr& c = ca[constraints[i]];
    if (c.isMarkedForDelete()) {
      c.freeUp();  // free up indirectly owned memory before implicitly deleting c during garbage collect
    } else {
      constraints[j++] = constraints[i];
    }
  }
  constraints.resize(j);
  if ((double)ca.wasted / (double)ca.at > 0.2) garbage_collect();
}

long Solver::clearLearnedConstraints() {
  assert(decisionLevel() == 0);

  std::vector<CRef> learnts;
  learnts.reserve(constraints.size() / 2);

  for (CRef& cr : constraints) {
    Constr& C = ca[cr];
    if (C.isMarkedForDelete() || external.count(C.id)) continue;
    BigVal eval = -C.degree();
    for (int j = 0; j < (int)C.size() && eval < 0; ++j)
      if (isUnit(Level, C.lit(j))) eval += C.coef(j);
    if (eval >= 0)
      removeConstraint(C, true);  // remove constraints satisfied at root
    else if (!C.isLocked()) {
      if (C.size() > 2 && C.lbd() > 2) learnts.push_back(cr);  // Keep all binary clauses and short LBDs
    }
  }

  for (size_t i = 0; i < learnts.size(); ++i) removeConstraint(ca[learnts[i]]);

  for (Lit l = -n; l <= n; ++l)
    for (int i = 0; i < (int)adj[l].size(); ++i) {
      if (ca[adj[l][i].cref].isMarkedForDelete()) aux::swapErase(adj[l], i--);
    }

  size_t i = 0;
  size_t j = 0;
  for (; i < constraints.size(); ++i) {
    Constr& c = ca[constraints[i]];
    if (c.isMarkedForDelete()) {
      c.freeUp();  // free up indirectly owned memory before implicitly deleting c during garbage collect
    } else {
      constraints[j++] = constraints[i];
    }
  }
  constraints.resize(j);
  if ((double)ca.wasted / (double)ca.at > 0.2) garbage_collect();
  return learnts.size();
}

// ---------------------------------------------------------------------
// Solving

double Solver::luby(double y, int i) const {
  // Find the finite subsequence that contains index 'i', and the
  // size of that subsequence:
  int size, seq;
  for (size = 1, seq = 0; size < i + 1; seq++, size = 2 * size + 1) {
  }
  while (size != i + 1) {
    size = (size - 1) >> 1;
    --seq;
    assert(size != 0);
    i = i % size;
  }
  return std::pow(y, seq);
}

bool Solver::checkSAT() {
  for (CRef cr : constraints) {
    const Constr& C = ca[cr];
    if (C.getOrigin() == Origin::FORMULA && !C.toExpanded(cePools)->isSatisfied(getLevel())) return false;
  }
  return true;
}

Lit Solver::pickBranchLit() {
  Var next = 0;
  // Activity based decision:
  while (next == 0 || !isUnknown(Pos, next)) {
    if (order_heap.empty())
      return 0;
    else
      next = order_heap.removeMax();
  }
  assert(phase[0] == 0);
  assert(lastSol[0] == 0);
  return phase[next];
}

void Solver::presolve() {
  firstRun = false;
  if (lpSolver) aux::timeCall<void>([&] { lpSolver->inProcess(); }, stats.LPTOTALTIME);
}

int Solver::varLevel(const int v) {
    if (Level[v] >= INF || Level[-v] >= INF)
      return std::min(Level[v], Level[-v]);
    return std::max(Level[v], Level[-v]);
}

bool Solver::varDefined(const int v) {
    return (Level[v] != INF or Level[-v] != INF);
}

bool Solver::varDecision(const int v) {
    return (Reason[v] == CRef_Undef);// and Reason[-v] == CRef_Undef);
}

inline void Solver::af1_mark_conflict(CRef confl, Var skipVar, int &n_toR, std::vector<int>& out_conflict, std::unordered_map<int, char>& seen) {
    //helper for analyzeFinal1--marks literals in conflict, adding some to out_conflict
    auto& c = ca[confl];
    for (int i = 0; i < c.size(); i++) {
        Lit l = c.lit(i); Var y = toVar(l);
        if(y == skipVar || varLevel(y) == 0 || seen.count(y)) //level zero lits not needed, skip already processed
            continue;
        if(Reason[y] == CRef_Undef) {
            if (Level[l] != INF)
              out_conflict.push_back(-l);
            else if (Level[-l] != INF)
              out_conflict.push_back(l);
            seen[y] = 2;
        }
        else { //must resolve later
            seen[y] = 1;
            ++n_toR;
        }
    }
}

inline void Solver::af1_mark_conflict(CeSuper confl, int &n_toR, std::vector<int>& out_conflict, std::unordered_map<int, char>& seen) {
    //helper for analyzeFinal1--marks literals in conflict, adding some to out_conflict
    for (Var y : confl->vars) {
        Lit l = confl->getLit(y);
        if(varLevel(y) == 0 || seen.count(y)) //level zero lits not needed, skip already processed
            continue;
        if(Reason[y] == CRef_Undef) {
            if (Level[l] != INF)
              out_conflict.push_back(-l);
            else if (Level[-l] != INF)
              out_conflict.push_back(l);
            seen[y] = 2;
        }
        else { //must resolve later
            seen[y] = 1;
            ++n_toR;
        }
    }
}

void Solver::analyzeFinal(CeSuper confl, std::vector<int> &out_core) {

  out_core.clear();

  if (decisionLevel() == 0) return;

  int n_toR {0}; //keep track of number of ltierals still to resolve.
  assert(decisionLevel() >= 1);
  std::unordered_map<int, char> seen;

  //uses seen markers to store the literals that need to be resolved
  //out of the conflict. Also uses seen markers to ensure that duplicate
  //lits not added to out_conflict.

  af1_mark_conflict(confl, n_toR, out_core, seen);
  //now do resolutions required to obtain final conflict.
  for (int i = trail.size()-1;  n_toR > 0 && i >= trail_lim[0]; i--) {
      Var x = toVar(trail[i]);
      if (Reason[x] == CRef_Undef) {
        if (seen.count(x) && seen[x] == 1)
          seen.erase(x);
        continue;
      }
      if(seen.count(x)) {
          if(seen[x] == 1) {
              //x has a reason and to be resolved away.
              --n_toR;
              af1_mark_conflict(Reason[x], x, n_toR, out_core, seen);
          }
          seen.erase(x);
      }
  }
}

void Solver::analyzeFinal(int p, std::vector<int> &out_core) {

  out_core.clear();
  out_core.push_back(-p);

  if (decisionLevel() == 0) return;

  std::unordered_set<int> seen;

  seen.insert(abs(p));

  for (int i = trail.size() - 1; i >= trail_lim[0]; i--) {
    int x = abs(trail[i]);
    if (seen.count(x)) {
      //std::cout << x << " was seen" << std::endl;
      if (varDecision(x)) {
        //std::cout << x << " decision " << std::endl;
        assert(varLevel(abs(trail[i])) < INF); 
        out_core.push_back(-trail[i]);
      } else {
        auto &c = ca[Reason[x]];
        for (int j = 0; j < (int)c.size(); j++) {
          //std::cout << c.lit(j) <<"@"<< varLevel(abs(c.lit(j))) << " ";
          if (varDefined(abs(c.lit(j)))) {
            if (c.coef(j) > 0 and Level[-c.lit(j)])
              seen.insert(abs(c.lit(j)));
            else if (c.coef(j) < 0 and Level[c.lit(j)])
              seen.insert(abs(c.lit(j)));
            else {
              //std::cout << "core reduced?!" << std::endl; // assignment does not reduce slack
            }
          } else { // else by prop on level 0 or unassigned
            //std::cout << "core reduced!" << std::endl; 
          }  
        }
        //std::cout << std::endl;
        seen.erase(x);
      }
    }
  }
}

thread_local std::random_device rd;
thread_local std::mt19937 g(rd());
// TODO: use a coroutine / yield instead of a SolveState return value
SolveAnswer Solver::solve(std::vector<Lit> assumptions, unsigned n_shuffles) {
  backjumpTo(0);
  bool runLP = firstRun || (lastState != SolveState::RESTARTED && lastState != SolveState::INPROCESSED);
  if (firstRun) presolve();

  while (true) {
    if (asynch_interrupt) throw asynchInterrupt;
    CeSuper confl = aux::timeCall<CeSuper>([&] { return runPropagation(runLP); }, stats.PROPTIME);
    runLP = !confl;
    if (confl) {
      assert(confl->hasNegativeSlack(Level));
      vDecayActivity();
      cDecayActivity();
      stats.NCONFL++;
      nconfl_to_restart--;
      if (stats.NCONFL % 1000 == 0 && options.verbosity.get() > 0) {
        printf("c #Conflicts: %10lld | #Constraints: %10lld\n", stats.NCONFL, (long long)constraints.size());
        if (options.verbosity.get() > 2) {
          // memory usage
          std::cout << "c total constraint space: " << ca.cap * 4 / 1024. / 1024. << "MB" << std::endl;
          std::cout << "c total #watches: ";
          long long cnt = 0;
          for (Lit l = -n; l <= n; l++) cnt += (long long)adj[l].size();
          std::cout << cnt << std::endl;
        }
      }
      if (decisionLevel() == 0) {
        if (logger) {
          confl->logInconsistency(Level, Pos, stats);
        }
        return {SolveState::UNSAT, {}, lastSol};
      }
      CeSuper analyzed = aux::timeCall<CeSuper>([&] { return analyze(confl); }, stats.CATIME);
      assert(analyzed);
      assert(analyzed->hasNegativeSlack(getLevel()));
      assert(analyzed->isSaturated());
      if (learnedStack.size() > 0 && learnedStack.back()->orig == Origin::FARKAS)
        learnConstraint(analyzed, Origin::LEARNEDFARKAS);  // TODO: ugly hack
      else
        learnConstraint(analyzed, Origin::LEARNED);
      
      /*if (decisionLevel() < (int)assumptions.size() && assumptions.size() > 0) {               //Conflict among assumptions
        std::vector<int> core;
        analyzeFinal(confl, core);
        backjumpTo(0);

        cores.push_back(core);

        if (shuffle_i >= n_shuffles) {
          return {SolveState::UNSAT, {}, lastSol};
        }

        shuffle_i++;
        std::random_shuffle(assumptions.begin(), assumptions.end());
      }*/
    } else {  // no conflict
      if (nconfl_to_restart <= 0) {
        backjumpTo(0);
        double rest_base = luby(options.lubyBase.get(), ++stats.NRESTARTS);
        nconfl_to_restart = (long long)rest_base * options.lubyMult.get();
        //        return {SolveState::RESTARTED, {}, lastSol}; // avoid this overhead for now
      }
      if (stats.NCONFL >= (stats.NCLEANUP + 1) * nconfl_to_reduce) {
        ++stats.NCLEANUP;
        if (options.verbosity.get() > 0) puts("c INPROCESSING");
        reduceDB();
        while (stats.NCONFL >= stats.NCLEANUP * nconfl_to_reduce) nconfl_to_reduce += options.dbCleanInc.get();
        if (lpSolver) aux::timeCall<void>([&] { lpSolver->inProcess(); }, stats.LPTOTALTIME);
        lastState = SolveState::INPROCESSED;
        return {SolveState::INPROCESSED, {}, lastSol};
      }
      
      Lit next = 0;
      if (decisionLevel() < (int)assumptions.size()) {
        next = assumptions[decisionLevel()];
        if (Level[-next] != INF) {
          std::vector<int> core;
          analyzeFinal(next, core);
          backjumpTo(0);

          cores.push_back(core);

          if (shuffle_i >= n_shuffles) {
            return {SolveState::UNSAT, {}, lastSol};
          }

          shuffle_i++;
          std::shuffle(assumptions.begin(), assumptions.end(), g);
        } 
        else if (Level[next] != INF) {
          trail_lim.push_back(trail.size());
          trail.push_back(0);//          Minisat doesn't have this
        } else {
          trail_lim.push_back(trail.size());
          uncheckedEnqueue(next, CRef_Undef);
        }
      } else {
        ++stats.NDECIDE;
        next = pickBranchLit();
        if (next == 0) {
          assert(order_heap.empty());
          assert((int)trail.size() == getNbVars());
          assert(checkSAT());
          lastSol.clear();
          lastSol.resize(getNbVars() + 1);
          lastSol[0] = 0;
          for (Var v = 1; v <= getNbVars(); ++v) lastSol[v] = isTrue(Level, v) ? v : -v;
          backjumpTo(0);
          return {SolveState::SAT, {}, lastSol};
        }
        decide(next);
        //cout << "decide " << next << " at " << decisionLevel() << endl;
      }
    }
  }
}

}  // namespace rs
