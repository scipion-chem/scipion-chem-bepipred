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
import os, subprocess, json

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
		cls._defineEmVar(BEPIPRED_DIC['home'], cls.getBepiPredDir())

	@classmethod
	def defineBinaries(cls, env):
		# todo: comprobar si esta ya instalado (requerimos definidos en el conf el comando de activacion, sea conda o python, y el home donde estan los scripts)
		# todo: si no hay comando de activacion, pero si home, comprobamos si es la carpeta y si es asi hacemos la instalacion.
		# todo: si no hay home definido ni se encuentra en el software em, no se instala nada (fallaria), ya se dirá que está mal instalado en el validateInstallation
		"""This function defines the binaries for each package."""
		if os.path.exists(cls.getVar(BEPIPRED_DIC['home'])):
			if not cls.getVar(BEPIPRED_DIC['activation']):
				cls._addBepiPredPackage(env, bepiHome=cls.getVar(BEPIPRED_DIC['home']))
			else:
				print('Environment activation command and HOME variable for BepiPred already found, instalaltion no needed')
		elif os.path.exists(cls.getVar(BEPIPRED_DIC['zip'])):
			cls._addBepiPredPackage(env, zipPath=cls.getVar(BEPIPRED_DIC['zip']))
		else:
			print()


	@classmethod
	def _addBepiPredPackage(cls, env, bepiHome=None, zipPath=None, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		BEPIPRED_INSTALLED = '%s_installed' % BEPIPRED_DIC['name']

		installationCmd = ''
		if not bepiHome and zipPath:
			bepiHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(BEPIPRED_DIC))
			installationCmd += f'unzip {zipPath} -d {bepiHome} && '

		installationCmd += f"cd {bepiHome} && sed -i 's/^torch==/#torch==/g' requirements.txt && "
		installationCmd += f"conda create -y -n {cls.getEnvName(BEPIPRED_DIC)} " \
											 f"python=3.9 --file requirements.txt && "
		installationCmd += f"touch {BEPIPRED_INSTALLED}"

		env.addPackage(BEPIPRED_DIC['name'], version=BEPIPRED_DIC['version'], tar=os.path.split(bepiHome)[-1],
									 commands=[(installationCmd, os.path.join(bepiHome, BEPIPRED_INSTALLED))],
									 neededProgs=["conda"], default=default, buildDir=bepiHome)

		cls._defineEmVar(BEPIPRED_DIC['activation'], cls.getEnvActivationCommand(BEPIPRED_DIC))


	@classmethod
	def validateInstallation(cls):
		""" Check if the installation of this protocol is correct.
         Returning an empty list means that the installation is correct
         and there are not errors. If some errors are found, a list with
         the error messages will be returned.
         """
		missingPaths = []
		if not os.path.exists(os.path.expanduser(cls.getVar(BEPIPRED_DIC['home']))):
			missingPaths.append("Path of BepiPred does not exist (%s) : %s " % (BEPIPRED_DIC['home'],
																																				 cls.getVar(BEPIPRED_DIC['home'])))
		return missingPaths

	@classmethod
	def getBepiPredDir(cls, fn=""):
		emDir = emConfig.EM_ROOT
		for file in os.listdir(emDir):
			if 'BepiPred' in file:
				return os.path.join(emDir, file, fn)
		print(f'BepiPred software could not be found in SOFTWARE directory ({emDir})')

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

