# pbo-ls-with-oracle

The solver is forked from [here](https://bitbucket.org/coreo-group/pbo-ihs-solver/src/master/).

## Dependencies
- Python version >= 3.6.9
- pybind11 version >= 2.5.0
- GCC version >= 8.4.0
- CPLEX version >= 12.8
- C++17 (i.e., a reasonably recent compiler)
- Boost library: https://www.boost.org

## Set up Python environment
- Create a new Python virtual environment with `python3 -m venv python-venv`
- Activate the environment with `source ./python-venv/bin/activate`
- Upgrade pip with `pip3 install --upgrade pip`
- Install requirements with `pip3 install -r ./requirements.txt`

## Building roundingsat solver
- Make sure to first activate the Python venv
- cd solver/roundingsat
- Follow instructions in solver/roundingsat/README.md file, chapter Compilation
- Alternatively follow instructions in chapter SoPlex if you want to utilise the LP relaxation feature
- mv roundingsat.cpython-36m-x86_64-linux-gnu.so ..
- Notice that roundingsat library file may have a different suffix

## Running the solver
- Run the solver with `python3 ./solver/solver.py path_to_instance.opb`