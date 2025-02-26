#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import sys
import os
from collections import OrderedDict

try:
    from distutils.core import setup
    from distutils.extension import Extension
except ImportError:
    print('The distutils package is required to build or install RMG Py.')
    raise
    
try:
    from Cython.Build import cythonize
    from Cython.Compiler import Options
except ImportError:
    print('Cython (http://www.cython.org/) is required to build or install RMG Py.')
    raise
    
try:
    import numpy
except ImportError:
    print('NumPy (http://numpy.scipy.org/) is required to build or install RMG Py.')
    raise

# Create annotated HTML files for each of the Cython modules
Options.annotate = True

directives = {
    # Set input language version to python 3
    'language_level': 3,
    # Turn on profiling capacity for all Cython modules
    # 'profile': True,
    # Embed call signatures in cythonized files - enable when building documentation
    # 'embedsignature': True,
}

################################################################################

main_ext_modules = [
    # RMG
    Extension('rmgpy.rmgobject', ['rmgpy/rmgobject.pyx']),
    # Kinetics
    Extension('rmgpy.kinetics.arrhenius', ['rmgpy/kinetics/arrhenius.pyx']),
    Extension('rmgpy.kinetics.chebyshev', ['rmgpy/kinetics/chebyshev.pyx']),
    Extension('rmgpy.kinetics.kineticsdata', ['rmgpy/kinetics/kineticsdata.pyx']),
    Extension('rmgpy.kinetics.falloff', ['rmgpy/kinetics/falloff.pyx']),
    Extension('rmgpy.kinetics.model', ['rmgpy/kinetics/model.pyx']),
    Extension('rmgpy.kinetics.tunneling', ['rmgpy/kinetics/tunneling.pyx']),
    Extension('rmgpy.kinetics.surface', ['rmgpy/kinetics/surface.pyx']),
    Extension('rmgpy.kinetics.uncertainties', ['rmgpy/kinetics/uncertainties.pyx']),
    # Molecules and molecular representations
    Extension('rmgpy.molecule.atomtype', ['rmgpy/molecule/atomtype.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.element', ['rmgpy/molecule/element.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.graph', ['rmgpy/molecule/graph.pyx'], include_dirs=['.']),
    Extension('rmgpy.molecule.group', ['rmgpy/molecule/group.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.molecule', ['rmgpy/molecule/molecule.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.symmetry', ['rmgpy/molecule/symmetry.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.vf2', ['rmgpy/molecule/vf2.pyx'], include_dirs=['.']),
    Extension('rmgpy.molecule.converter', ['rmgpy/molecule/converter.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.translator', ['rmgpy/molecule/translator.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.util', ['rmgpy/molecule/util.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.inchi', ['rmgpy/molecule/inchi.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.resonance', ['rmgpy/molecule/resonance.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.pathfinder', ['rmgpy/molecule/pathfinder.py'], include_dirs=['.']),
    Extension('rmgpy.molecule.kekulize', ['rmgpy/molecule/kekulize.pyx'], include_dirs=['.']),
    # Pressure dependence
    Extension('rmgpy.pdep.collision', ['rmgpy/pdep/collision.pyx']),
    Extension('rmgpy.pdep.configuration', ['rmgpy/pdep/configuration.pyx']),
    Extension('rmgpy.pdep.me', ['rmgpy/pdep/me.pyx']),
    Extension('rmgpy.pdep.msc', ['rmgpy/pdep/msc.pyx']),
    Extension('rmgpy.pdep.reaction', ['rmgpy/pdep/reaction.pyx']),
    Extension('rmgpy.pdep.rs', ['rmgpy/pdep/rs.pyx']),
    Extension('rmgpy.pdep.cse', ['rmgpy/pdep/cse.pyx']),
    # Statistical mechanics
    Extension('rmgpy.statmech.conformer', ['rmgpy/statmech/conformer.pyx']),
    Extension('rmgpy.statmech.mode', ['rmgpy/statmech/mode.pyx']),
    Extension('rmgpy.statmech.rotation', ['rmgpy/statmech/rotation.pyx']),
    Extension('rmgpy.statmech.schrodinger', ['rmgpy/statmech/schrodinger.pyx']),
    Extension('rmgpy.statmech.torsion', ['rmgpy/statmech/torsion.pyx']),
    Extension('rmgpy.statmech.translation', ['rmgpy/statmech/translation.pyx']),
    Extension('rmgpy.statmech.vibration', ['rmgpy/statmech/vibration.pyx']),
    # Thermodynamics
    Extension('rmgpy.thermo.thermodata', ['rmgpy/thermo/thermodata.pyx']),
    Extension('rmgpy.thermo.model', ['rmgpy/thermo/model.pyx']),
    Extension('rmgpy.thermo.nasa', ['rmgpy/thermo/nasa.pyx']),
    Extension('rmgpy.thermo.wilhoit', ['rmgpy/thermo/wilhoit.pyx']),
    # Miscellaneous
    Extension('rmgpy.constants', ['rmgpy/constants.py'], include_dirs=['.']),
    Extension('rmgpy.quantity', ['rmgpy/quantity.py'], include_dirs=['.']),
    Extension('rmgpy.reaction', ['rmgpy/reaction.py'], include_dirs=['.']),
    Extension('rmgpy.species', ['rmgpy/species.py'], include_dirs=['.']),
    Extension('rmgpy.chemkin', ['rmgpy/chemkin.pyx'], include_dirs=['.']),
]

solver_ext_modules = [
    Extension('rmgpy.solver.base', ['rmgpy/solver/base.pyx'], include_dirs=['.']),
    Extension('rmgpy.solver.simple', ['rmgpy/solver/simple.pyx'], include_dirs=['.']),
    Extension('rmgpy.solver.liquid', ['rmgpy/solver/liquid.pyx'], include_dirs=['.']),
    Extension('rmgpy.solver.mbSampled', ['rmgpy/solver/mbSampled.pyx'], include_dirs=['.']),
    Extension('rmgpy.solver.surface', ['rmgpy/solver/surface.pyx'], include_dirs=['.']),
]

arkane_ext_modules = [
    # RMG
    Extension('rmgpy.rmgobject', ['rmgpy/rmgobject.pyx']),
    # Kinetics
    Extension('rmgpy.kinetics.arrhenius', ['rmgpy/kinetics/arrhenius.pyx']),
    Extension('rmgpy.kinetics.chebyshev', ['rmgpy/kinetics/chebyshev.pyx']),
    Extension('rmgpy.kinetics.kineticsdata', ['rmgpy/kinetics/kineticsdata.pyx']),
    Extension('rmgpy.kinetics.falloff', ['rmgpy/kinetics/falloff.pyx']),
    Extension('rmgpy.kinetics.model', ['rmgpy/kinetics/model.pyx']),
    Extension('rmgpy.kinetics.tunneling', ['rmgpy/kinetics/tunneling.pyx']),
    # Pressure dependence
    Extension('rmgpy.pdep.collision', ['rmgpy/pdep/collision.pyx']),
    Extension('rmgpy.pdep.configuration', ['rmgpy/pdep/configuration.pyx']),
    Extension('rmgpy.pdep.me', ['rmgpy/pdep/me.pyx']),
    Extension('rmgpy.pdep.msc', ['rmgpy/pdep/msc.pyx']),
    Extension('rmgpy.pdep.reaction', ['rmgpy/pdep/reaction.pyx']),
    Extension('rmgpy.pdep.rs', ['rmgpy/pdep/rs.pyx']),
    Extension('rmgpy.pdep.cse', ['rmgpy/pdep/cse.pyx']),
    # Statistical mechanics
    Extension('rmgpy.statmech.conformer', ['rmgpy/statmech/conformer.pyx']),
    Extension('rmgpy.statmech.mode', ['rmgpy/statmech/mode.pyx']),
    Extension('rmgpy.statmech.rotation', ['rmgpy/statmech/rotation.pyx']),
    Extension('rmgpy.statmech.schrodinger', ['rmgpy/statmech/schrodinger.pyx']),
    Extension('rmgpy.statmech.torsion', ['rmgpy/statmech/torsion.pyx']),
    Extension('rmgpy.statmech.translation', ['rmgpy/statmech/translation.pyx']),
    Extension('rmgpy.statmech.vibration', ['rmgpy/statmech/vibration.pyx']),
    # Thermodynamics
    Extension('rmgpy.thermo.thermodata', ['rmgpy/thermo/thermodata.pyx']),
    Extension('rmgpy.thermo.model', ['rmgpy/thermo/model.pyx']),
    Extension('rmgpy.thermo.nasa', ['rmgpy/thermo/nasa.pyx']),
    Extension('rmgpy.thermo.wilhoit', ['rmgpy/thermo/wilhoit.pyx']),
    # Miscellaneous
    Extension('rmgpy.constants', ['rmgpy/constants.py'], include_dirs=['.']),
    Extension('rmgpy.quantity', ['rmgpy/quantity.py'], include_dirs=['.']),
]

################################################################################

ext_modules = []
if 'install' in sys.argv:
    # This is so users can still do simply `python setup.py install`
    ext_modules.extend(main_ext_modules)
    ext_modules.extend(solver_ext_modules)
if 'main' in sys.argv:
    # This is for `python setup.py build_ext main`
    sys.argv.remove('main')
    ext_modules.extend(main_ext_modules)
if 'solver' in sys.argv:
    # This is for `python setup.py build_ext solver`
    sys.argv.remove('solver')
    ext_modules.extend(solver_ext_modules)
if 'arkane' in sys.argv:
    # This is for `python setup.py build_ext arkane`
    sys.argv.remove('arkane')
    ext_modules.extend(main_ext_modules)
    ext_modules.extend(arkane_ext_modules)
if 'minimal' in sys.argv:
    # This starts with the full install list, but removes anything that has a pure python mode
    # i.e. in only includes things whose source is .pyx
    sys.argv.remove('minimal')
    temporary_list = []
    temporary_list.extend(main_ext_modules)
    temporary_list.extend(solver_ext_modules)
    for module in temporary_list:
        for source in module.sources:
            if os.path.splitext(source)[1] == '.pyx':
                ext_modules.append(module)

# Remove duplicates while preserving order:
ext_modules = list(OrderedDict.fromkeys(ext_modules))

scripts = [
    'Arkane.py',
    'rmg.py',
    'scripts/checkModels.py',
    'scripts/convertFAME.py',
    'scripts/diffModels.py',
    'scripts/generateChemkinHTML.py',
    'scripts/generateFluxDiagram.py',
    'scripts/generateReactions.py',
    'scripts/machineWriteDatabase.py',
    'scripts/mergeModels.py',
    'scripts/simulate.py',
    'scripts/standardizeModelSpeciesNames.py',
    'scripts/thermoEstimator.py',
    'testing/databaseTest.py',
]

modules = []
for root, dirs, files in os.walk('rmgpy'):
    if 'test_data' in root:
        continue
    for f in files:
        if f.endswith('.py') or f.endswith('.pyx'):
            if 'Test' not in f and '__init__' not in f:
                module = 'rmgpy' + root.partition('rmgpy')[-1].replace('/', '.') + '.' + f.partition('.py')[0]
                modules.append(module)
for root, dirs, files in os.walk('arkane'):
    if 'data' in root:
        continue
    for f in files:
        if f.endswith('.py') or f.endswith('.pyx'):
            if 'Test' not in f and '__init__' not in f:
                module = 'arkane' + root.partition('arkane')[-1].replace('/', '.') + '.' + f.partition('.py')[0]
                modules.append(module)

# Read the version number
exec(open('rmgpy/version.py').read())

# Initiate the build and/or installation
setup(
    name='RMG-Py',
    version=__version__,
    description='Reaction Mechanism Generator',
    author='William H. Green and the RMG Team',
    author_email='rmg_dev@mit.edu',
    url='http://reactionmechanismgenerator.github.io',
    packages=['rmgpy', 'arkane'],
    py_modules=modules,
    scripts=scripts,
    ext_modules=cythonize(ext_modules, build_dir='build', compiler_directives=directives),
    include_dirs=['.', numpy.get_include()],
)
