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

#include <csignal>
#include <fstream>
#include <memory>
#include "aux.hpp"
#include "globals.hpp"
#include "parsing.hpp"
#include "run.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace rs {

bool asynch_interrupt;
Options options;
Stats stats;

}  // namespace rs

static void SIGINT_interrupt([[maybe_unused]] int signum) { rs::asynch_interrupt = true; }

static void SIGINT_exit([[maybe_unused]] int signum) {
  printf("\n*** INTERRUPTED ***\n");
  exit(1);
}

class Roundingsat {
  public:
  void init(int n) {
    rs::run::solver.init();
    rs::run::solver.setNbVars(std::max(n * 2, 5000), true);
  }

  void addConstraint(py::dict expression, long rhs) {
    rs::CeArb input = rs::run::solver.cePools.takeArb();
    input->reset();
    for (auto pair : expression) {
      input->addLhs(pair.second.cast<int>(), pair.first.cast<int>());
    }
    input->addRhs(rhs);
    if (rs::run::solver.addConstraint(input, rs::Origin::FORMULA).second == rs::ID_Unsat) rs::quit::exit_UNSAT(rs::run::solver);
  }

  //int main(int argc, char** argv) {
  void solve(py::list assums, unsigned n_shuffles) {
    rs::stats.STARTTIME = rs::aux::cpuTime();
    rs::asynch_interrupt = false;

    signal(SIGINT, SIGINT_exit);
    signal(SIGTERM, SIGINT_exit);
    signal(SIGXCPU, SIGINT_exit);
    signal(SIGINT, SIGINT_interrupt);
    signal(SIGTERM, SIGINT_interrupt);
    signal(SIGXCPU, SIGINT_interrupt);

    //rs::options.parseCommandLine(argc, argv);

    //if (rs::options.verbosity.get() > 0) {
    //  std::cout << "c RoundingSat 2\n";
    //  std::cout << "c branch " << EXPANDED(GIT_BRANCH) << "\n";
    //  std::cout << "c commit " << EXPANDED(GIT_COMMIT_HASH) << std::endl;
    //}

    //if (!instance_file.empty()) {
    //  std::ifstream fin(instance_file);
    //  if (!fin) rs::quit::exit_ERROR({"Could not open ", instance_file});
    //  rs::parsing::file_read(fin, rs::run::solver, objective);
    /*} else {
      if (rs::options.verbosity.get() > 0) std::cout << "c No filename given, reading from standard input" << std::endl;
      rs::parsing::file_read(std::cin, rs::run::solver, objective);
    }*/

    rs::CeArb objective = rs::run::solver.cePools.takeArb();
    rs::run::solver.initLP(objective);

    std::vector<int> assumptions;
    for (auto l : assums) {
      assumptions.push_back(l.cast<int>());
    }
    rs::run::solver.shuffle_i = 1;
    rs::run::solver.cores.clear();

    rs::run::run(objective, assumptions, n_shuffles);
  }

  py::tuple getResult() {
    py::list solution;
    for (auto s : rs::run::solver.lastSol)
      if (s != 0)
        solution.append(s);
    py::list cores;
    for (const std::vector<int>& c : rs::run::solver.cores) {
      py::list core;
      for (auto cc : c)
        core.append(cc);
      cores.append(core);
    }  
    return py::make_tuple(cores, solution, (int)rs::run::solver.lastState + 1);
  }

  long clearLearnedConstraints() {
    return rs::run::solver.clearLearnedConstraints();
  }

  void print() {
    rs::stats.print();
  }
};

PYBIND11_MODULE(roundingsat, m) {
    py::class_<Roundingsat>(m, "Roundingsat")
        .def(py::init<>())
        .def("init", &Roundingsat::init)
        //.def_readwrite("verbosity", &Roundingsat::verbosity)
        .def("getResult", &Roundingsat::getResult)
        .def("addConstraint", &Roundingsat::addConstraint)
        .def("solve", &Roundingsat::solve)
        .def("clearLearnedConstraints", &Roundingsat::clearLearnedConstraints)
        .def("print", &Roundingsat::print);
}
