# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using ConPLex models for virtual screening
"""

# General imports
import os, subprocess, json, pathlib

# Scipion em imports
import pwem
from pwem import Config as emConfig
from scipion.install.funcs import InstallHelper

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'dtu_logo.jpeg'

class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineVar(BEPIPRED_DIC['home'], cls.getBepiPredDir())
		cls._defineVar(BEPIPRED_DIC['activation'], cls.getEnvActivationCommand(BEPIPRED_DIC))
		cls._defineVar(BEPIPRED_DIC['zip'], None)

	@classmethod
	def defineBinaries(cls, env):
		"""This function defines the binaries for each package."""
		if cls.checkVarPath('zip'):
			cls._addBepiPredPackage(env, zipPath=cls.getVar(BEPIPRED_DIC['zip']))
		elif cls.checkVarPath('home'):
			if not cls.checkCallEnv(BEPIPRED_DIC):
				cls._addBepiPredPackage(env, bepiHome=cls.getVar(BEPIPRED_DIC['home']))
			# else:
			# 	print('Environment activation command and HOME variables for BepiPred already found, installation no needed')

	@classmethod
	def _addBepiPredPackage(cls, env, bepiHome=None, zipPath=None, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		BEPIPRED_INSTALLED = '%s_installed' % BEPIPRED_DIC['name']

		installationCmd = ''
		if not bepiHome and zipPath:
			bepiHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(BEPIPRED_DIC))
			installationCmd += f'unzip -q {zipPath} -d {bepiHome} && ' \
												 f'mv {bepiHome}/BepiPred3_src/* {bepiHome} && rm -r {bepiHome}/BepiPred3_src && '

		installationCmd += f"cd {bepiHome} && sed -i 's/^torch==/#torch==/g' requirements.txt && "
		installationCmd += f"conda create -y -n {cls.getEnvName(BEPIPRED_DIC)} " \
											 f"python=3.9 --file requirements.txt && "
		installationCmd += f"touch {BEPIPRED_INSTALLED}"

		env.addPackage(BEPIPRED_DIC['name'], version=BEPIPRED_DIC['version'],
									 commands=[(installationCmd, os.path.join(bepiHome, BEPIPRED_INSTALLED))], tar='void.tgz',
									 neededProgs=["conda"], default=default, buildDir=os.path.split(bepiHome)[-1])


	@classmethod
	def validateInstallation(cls):
		""" Check if the installation of this protocol is correct. Returning an empty list means that the installation
		is correct and there are not errors. If some errors are found, a list with the error messages will be returned."""
		mPaths = []
		if not cls.checkVarPath('home'):
			mPaths.append(f"Path of BepiPred home (folder like BepiPred3_src) does not exist.\n"
										f"You must either define it in the scipion.conf (as {BEPIPRED_DIC['home']} = <pathToBepiPred_src>) "
										f"or define the location of the raw dowloaded ZIP file (like bepipred-3.0b.src.zip) as "
										f"{BEPIPRED_DIC['zip']} = <pathToBepiPredZip>.\nAlternatively, you can move the home folder into "
										f"{emConfig.EM_ROOT} keeping the '{bepiPattern}' pattern.")

		if not cls.checkCallEnv(BEPIPRED_DIC):
			mPaths.append(f"Activation of the BepiPred environment failed.\n")

		if len(mPaths) > 0:
			mPaths.append(NOINSTALL_WARNING)
		return mPaths

	@classmethod
	def getBepiPredDir(cls, fn=""):
		emDir = emConfig.EM_ROOT
		for file in os.listdir(emDir):
			if bepiPattern in file.lower():
				return os.path.join(emDir, file, fn)
		# print(f'BepiPred software could not be found in SOFTWARE directory ({emDir})')
		return os.path.join(emConfig.EM_ROOT, cls.getEnvName(BEPIPRED_DIC))

	@classmethod
	def checkVarPath(cls, var='home'):
		'''Check if a plugin variable exists and so do its path'''
		exists = False
		varValue = cls.getVar(BEPIPRED_DIC[var])
		if varValue and os.path.exists(varValue):
			exists = True
		return exists

	@classmethod
	def checkCallEnv(cls, packageDic):
		actCommand = cls.getVar(packageDic['activation'])
		try:
			if 'conda' in actCommand and not 'shell.bash hook' in actCommand:
				actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
			subprocess.check_output(actCommand, shell=True)
			envFine = True
		except subprocess.CalledProcessError as e:
			envFine = False
		return envFine


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runScript(cls, protocol, scriptName, args, envDict, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDict), 'python', scriptName)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	# ---------------------------------- Utils functions-----------------------

