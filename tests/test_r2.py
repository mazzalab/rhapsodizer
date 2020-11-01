""""
File description
"""

# import os

__author__ = "Tommaso Mazza"
__copyright__ = "Copyright 2020, The rhapsodizer Project"
__version__ = "0.0.1"
__maintainer__ = "Tommaso Mazza"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = "01/11/2020"
__creator__ = "t.mazza"
__license__ = u"""
  Copyright (C) 2020  Tommaso Mazza <t.mazza@css-mendel.it>
  Viale Regina Margherita 261, 00198 Rome, Italy

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
  """

from unittest import TestCase, TestLoader, TextTestRunner


class TestR2(TestCase):
    def test_read_st(self):
        self.fail()

    def test_read_cartridge_index(self):
        self.fail()

    def test_get_file_name(self):
        self.fail()

    def test_parse_bam(self):
        self.fail()

    def test_parse_fastq(self):
        self.fail()


suite = TestLoader().loadTestsFromTestCase(TestR2)
TextTestRunner(verbosity=2).run(suite)
