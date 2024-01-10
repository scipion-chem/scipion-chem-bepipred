# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, glob, shutil, subprocess

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.objects import SequenceChem, SetOfSequencesChem

from .. import Plugin as bepiPlugin
from ..constants import BEPIPRED_DIC

class ProtBepiPredPrediction(EMProtocol):
  """Run a prediction using BepiPred to extract B-cell epitopes"""
  _label = 'bepipred prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequences",
                    label='Input protein sequences: ',
                    help="Set of protein sequences to perform the screening on")

  def _insertAllSteps(self):
    self._insertFunctionStep(self.trialStep)

  def trialStep(self):
    for par in ['home', 'activation']:
      print(f'{BEPIPRED_DIC[par]}: {bepiPlugin.getVar(BEPIPRED_DIC[par])}')

