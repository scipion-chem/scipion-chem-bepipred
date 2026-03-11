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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs

from .. import Plugin as bepiPlugin
from ..constants import BEPIPRED_DIC

class ProtBepiPredPrediction(EMProtocol):
  """Run a prediction using BepiPred to extract B-cell epitopes

  AI Generated:
        ProtBepiPredPrediction - User Manual

        Overview
        --------
        The ProtBepiPredPrediction protocol predicts B-cell epitopes within a
        protein sequence using the BepiPred framework. B-cell epitopes are
        regions of an antigen that are recognized by antibodies and therefore
        play a crucial role in immune response, vaccine development, and
        immunological research.

        This protocol analyzes an input protein sequence and evaluates each
        residue according to its predicted probability of belonging to an
        epitope region. The predictions are generated using pretrained machine
        learning models that combine sequence-based features with protein
        embedding representations.

        Input Requirements
        ------------------
        - **Sequence**:
            A protein sequence provided as a `Sequence` object. The sequence
            should represent the antigen or protein of interest and must contain
            valid amino acid residues.

        Parameters
        ----------
        - **Prediction type (predType)**:
            Determines the BepiPred prediction model used:
                • `mjv_pred`: Majority vote ensemble model.
                • `vt_pred`: Variable threshold prediction based on ensemble
                  probabilities.

        - **Top proportion (top)**:
            Defines the proportion of residues with the highest prediction
            scores that will be considered epitope candidates.

        - **Add sequence length (addSeqLen)**:
            Adds sequence length information to the ESM-based encoding used by
            the prediction model.

        - **Extract linear epitopes (linearEp)**:
            Enables smoothing of prediction scores using a rolling window to
            identify contiguous epitope regions.

        - **Rolling window size (rWindow)**:
            Window size used to compute the rolling average of epitope
            probabilities when extracting linear epitopes.

        Epitope Extraction Parameters
        -----------------------------
        - **Average threshold (avThres)**:
            Minimum prediction score required for a residue to be considered
            epitope-positive.

        - **Soft threshold (useSoft)**:
            Allows inclusion of residues with slightly lower scores when they
            occur between positive residues within an epitope region.

        - **Soft threshold value (softThres)**:
            Defines the tolerance range for residues included by the soft
            threshold.

        - **Soft threshold size (nSoft)**:
            Maximum number of consecutive residues that can be included using
            the soft threshold rule.

        Epitope Size Constraints
        ------------------------
        - **Set epitope size limits (setSize)**:
            Enables filtering of predicted epitopes by length.

        - **Minimum epitope size (minSize)**:
            Minimum number of residues allowed in a predicted epitope.

        - **Maximum epitope size (maxSize)**:
            Maximum number of residues allowed in a predicted epitope.

        Workflow
        --------
        1. **Input preparation**
           - The protein sequence is exported to FASTA format.

        2. **BepiPred prediction**
           - The sequence is analyzed using the selected BepiPred model.
           - Residue-level epitope probability scores are generated.

        3. **Score processing**
           - Prediction scores are optionally smoothed using a rolling window.
           - Residues exceeding the threshold are grouped into candidate
             epitopes.

        4. **Epitope extraction**
           - Consecutive residues forming epitope segments are identified.
           - Soft threshold rules may include borderline residues.
           - Size constraints are applied if specified.

        5. **Output generation**
           - Each predicted epitope is converted into a `SequenceROI`
             representing a region of interest within the protein sequence.

        Outputs
        -------
        - **SetOfSequenceROIs**:
            A collection of predicted epitope regions associated with the input
            protein sequence. Each region contains:
                • Epitope sequence
                • Start and end positions in the protein
                • Annotation describing the predicted B-cell epitope

        - **Raw prediction files**:
            Intermediate CSV files containing residue-level prediction scores.

        Interpretation
        --------------
        - High prediction scores indicate residues more likely to belong to
          antibody-recognized regions.
        - Extracted epitopes represent contiguous segments with strong
          epitope probability signals.
        - Adjusting thresholds and smoothing parameters allows users to
          control sensitivity versus specificity of predictions.

        Practical Recommendations
        -------------------------
        - Use complete antigen sequences for best results.
        - Adjust threshold parameters depending on the desired number of
          predicted epitopes.
        - Enable linear epitope extraction when identifying contiguous
          antibody-binding segments.

        Warnings
        --------
        - Predictions are computational estimates and should be validated
          experimentally when possible.
        - Prediction accuracy may vary depending on protein structure and
          antigen context.
        - Extremely short or disordered sequences may produce less reliable
          results.

        Final Perspective
        -----------------
        ProtBepiPredPrediction provides an automated approach to identify
        potential B-cell epitopes directly from protein sequences. By combining
        machine learning predictions with flexible epitope extraction
        parameters, the protocol enables efficient screening of antigenic
        regions for vaccine design, antibody development, and immunological
        studies.

    """
  _label = 'bepipred prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                    label='Input protein sequence: ',
                    help="Protein sequence to perform the screening on")

    pGroup = form.addGroup('Parameters')
    pGroup.addParam('predType', params.EnumParam, label='Prediction type: ', default=0, choices=["mjv_pred", "vt_pred"],
                    help="Prediction model to use, either majority vote ensemble or variable threshold predicition "
                         "on average ensemble posistive probabilities")
    pGroup.addParam('top', params.FloatParam, label='Top proportion: ', default=0.2,
                    expertLevel=params.LEVEL_ADVANCED, help="Top proportion of epitope residues")
    pGroup.addParam('addSeqLen', params.BooleanParam, label='Add sequence length: ', default=False,
                    expertLevel=params.LEVEL_ADVANCED, help="Add sequence lengths to esm-encodings")

    pGroup.addParam('linearEp', params.BooleanParam, label='Extract linear epitopes: ', default=False,
                    help="Whether to pass a smoothing window over the results to extract linear epitopes")
    pGroup.addParam('rWindow', params.IntParam, label='Rolling window size: ', default=9, condition='linearEp',
                    help="Window size to use for rolling average on B-cell epitope probability scores.")

    eGroup = form.addGroup('Epitope extraction')
    eGroup.addParam('avThres', params.FloatParam, label='Average threshold: ', default=0.1512,
                    help="Threshold to use, when making predictions for considering a residue to be epitope positive")
    eGroup.addParam('useSoft', params.BooleanParam, label='Use a soft threshold: ', default=False,
                    expertLevel=params.LEVEL_ADVANCED,
                    help="Use a soft threshold when extracting the epitopes by score")
    eGroup.addParam('softThres', params.FloatParam, label='Soft threshold: ', default=0.1,
                    expertLevel=params.LEVEL_ADVANCED, condition='useSoft',
                    help="Soft threshold to use. If the score of a negative residue between two positive residues is "
                         "over the threshold-(threshold*softThreshold), it is also included as positive")
    eGroup.addParam('nSoft', params.IntParam, label='Soft threshold size: ', default=1,
                    expertLevel=params.LEVEL_ADVANCED, condition='useSoft',
                    help="Defines N as the number of residues for which the soft threshold can be applied sequentially."
                         "If the score of N negative residue between two positive residues are "
                         "over threshold-(threshold*softThreshold), they are also included as positive")

    eGroup.addParam('setSize', params.BooleanParam, label='Set epitope size limits: ', default=False,
                    help="Whether to establish minimum and maximum epitope size")
    eGroup.addParam('minSize', params.IntParam, label='Minimum epitope size: ', default=3, condition='setSize',
                    help="Minimum epitope size")
    eGroup.addParam('maxSize', params.IntParam, label='Maximum epitope size: ', default=15, condition='setSize',
                    help="Maximum epitope size")


  def _insertAllSteps(self):
    self._insertFunctionStep(self.bepipredStep)
    self._insertFunctionStep(self.createOutputStep)

  def writeInputFasta(self):
    faFile = self._getExtraPath('inputSequence.fa')
    self.inputSequence.get().exportToFile(faFile)
    return os.path.abspath(faFile)

  def bepipredStep(self):
    faFile = self.writeInputFasta()
    oDir = os.path.abspath(self._getExtraPath())

    bepiArgs = f'-i {faFile} -o {oDir} -esm_dir {bepiPlugin.getVar(BEPIPRED_DIC["home"])} ' \
               f'-pred {self.getEnumText("predType")} -t {self.avThres.get()} -top {self.top.get()} '
    if self.linearEp.get():
      bepiArgs += f'-rolling_window_size {self.rWindow.get()} '
    if self.addSeqLen.get():
      bepiArgs += '-add_seq_len '

    bepiPlugin.runBepiPred(self, bepiArgs)

  def createOutputStep(self):
    epiDic = self.parseResults(minLen=self.minSize.get(), maxLen=self.maxSize.get(), threshold=self.avThres.get(),
                               softThres=self.softThres.get(), softN=self.nSoft.get())

    inpSeq = self.inputSequence.get()
    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    for idxI, (epitope, score) in epiDic[list(epiDic.keys())[0]].items():
      idxs = [idxI, idxI+len(epitope)]
      roiSeq = Sequence(sequence=epitope, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                        description=f'BepiPred epitope')
      seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
      outROIs.append(seqROI)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def parseResults(self, minLen, maxLen, threshold=0.1512, softThres=0.1, softN=1):
    '''Parse the results in the raw_output.csv file generated by BepiPred and returns a dictionary
    as {protId: {position: [epitope, meanScore]}} with the epitopes passing the threshold for BepiPred-3.0 linear epitope score.
    To pass the threshold, an epitope needs to have a lenght between minLen and maxLen.
    Also, all residue scores must be over the threshold or, optionally, between these residues there can be softN residues
    that passes the soft threshold, defined as thres-thres*softThres. This allows sort some close to the threshold
    residues inside the epitope.
    '''
    oDir = self._getExtraPath()

    epiDic = {}
    curEp, iniEp = '', -1
    with open(os.path.join(oDir, 'raw_output.csv')) as f:
      f.readline()
      for line in f:
        protId, res, score3D, scoreLinear = line.split(',')
        score = scoreLinear if self.linearEp.get() else score3D
        if protId not in epiDic:
          epiDic[protId] = {}
          i = 1
          curEp, iniEp, curSoft, scores = '', -1, 0, []

        if float(score) >= threshold:
          # If the threshold is passed
          if not curEp:
            # Save protein position
            iniEp = i
          # Add residue to epitope, initilize soft threshold
          curEp += res
          scores.append(float(score))
          curSoft = 0

        elif self.useSoft.get() and curEp and float(score) >= (threshold - threshold * softThres) and curSoft < softN:
          # If we have started defining an epitope and the soft thershold is passed, add residue to epitope
          curSoft += 1
          curEp += res
          scores.append(float(score))

        else:
          # If nor threshold is passed
          if curSoft > 0:
            # If last added position(s) was soft threshold, remove them
            curEp = curEp[:-curSoft]
          if len(scores) > 0:
            if not self.setSize.get() or (len(curEp) >= minLen and len(curEp) <= maxLen):
              # Save current epitope if the len conditions are met
              epiDic[protId][iniEp] = [curEp, sum(scores)/len(scores)]
            # Initialize current epitope
            curEp, iniEp, curSoft, scores = '', -1, 0, []

        i += 1
    return epiDic
