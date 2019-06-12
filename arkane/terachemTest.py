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

import unittest
import os

from rmgpy.statmech import Conformer, IdealGasTranslation, LinearRotor, NonlinearRotor, HarmonicOscillator, \
    HinderedRotor

from arkane.terachem import TeraChemLog


################################################################################


class TeraChemTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.TeraChem module, used for reading
    and writing TeraChem files.
    """

    def testNumberOfAtomsFromTeraChemLog(self):
        """
        Uses a TeraChem log files to test that number of atoms can be properly read.
        """
        log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'terachem', 'minimize_output.out'))
        self.assertEqual(log.getNumberOfAtoms(), 6)
        log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'terachem', 'coord.xyz'))
        self.assertEqual(log.getNumberOfAtoms(), 6)
        log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'terachem', 'output.geometry'))
        self.assertEqual(log.getNumberOfAtoms(), 6)

    # def testEnergyFromTeraChemLog(self):
    #     """
    #     Uses a TeraChem log files to test that
    #     molecular energies can be properly read.
    #     """
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'npropyl.out'))
    #     self.assertAlmostEqual(log.loadEnergy(), -310896203.5432524, delta=1e-5)
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'co.out'))
    #     self.assertAlmostEqual(log.loadEnergy(), -297402545.0217114, delta=1e-5)
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'CH4_sp_TeraChem.out'))
    #     self.assertAlmostEqual(log.loadEnergy(), -106356735.53661588, delta=1e-5)
    #
    # def testLoadVibrationsFromTeraChemLog(self):
    #     """
    #     Uses a TeraChem log files to test that
    #     molecular energies can be properly read.
    #     """
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'npropyl.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #     self.assertEqual(len(conformer.modes[2]._frequencies.getValue()), 24)
    #     self.assertEqual(conformer.modes[2]._frequencies.getValue()[5], 881.79)
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'co.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #     self.assertEqual(len(conformer.modes[2]._frequencies.getValue()), 1)
    #     self.assertEqual(conformer.modes[2]._frequencies.getValue(), 2253.16)
    #
    # def testLoadNpropylModesFromTeraChemLog(self):
    #     """
    #     Uses a TeraChem log file for npropyl to test that its
    #     molecular modes can be properly read.
    #     """
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'npropyl.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)
    #
    # def testSpinMultiplicityFromTeraChemLog(self):
    #     """
    #     Uses a TeraChem log file for npropyl to test that its
    #     molecular degrees of freedom can be properly read.
    #     """
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'npropyl.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #     self.assertEqual(conformer.spinMultiplicity, 2)
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'co.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #     self.assertEqual(conformer.spinMultiplicity, 1)
    #
    # def testLoadCOModesFromTeraChemLog(self):
    #     """
    #     Uses a TeraChem log file for CO to test that its
    #     molecular degrees of freedom can be properly read.
    #     """
    #     log = TeraChemLog(os.path.join(os.path.dirname(__file__), 'data', 'co.out'))
    #     conformer, unscaled_frequencies = log.loadConformer()
    #     E0 = log.loadEnergy()
    #
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, IdealGasTranslation)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, LinearRotor)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, NonlinearRotor)]) == 0)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HarmonicOscillator)]) == 1)
    #     self.assertTrue(len([mode for mode in conformer.modes if isinstance(mode, HinderedRotor)]) == 0)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
