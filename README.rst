================================
BepiPred scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of tools provided by BepiPred.
This plugin allows to use programs from the BepiPred software
within the Scipion framework. **You need to download the BepiPred files
before installing the plugin, see section "Download BepiPred" for details**.

================================
Download BepiPred
================================

BepiPred is a software meant for academic use only. You can download BepiPred-3.0 in
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=bepipred&version=3.0&packageversion=3.0b&platform=src.

|

Once you obtain the software file (a zip) you have several options to help Scipion finding it:

Option 1) Edit the scipion.conf file and add the variable: BEPIPRED_ZIP = <PathToBepiPredZip>.
This way, Scipion will unzip and move the corresponding files to the scipion/software/em folder and install BepiPred.

Option 2) If you have unzipped BepiPred yourself you can either:

2.1) Move the folder (of the form BepiPred3_src) to the scipion/software/em folder. Scipion will find it there.

2.2) Specify the location of the BepiPred folder in the scipion.conf file as: BEPIPRED_HOME = <PathToBepiPred3_src>

Option 3) If you have already installed BepiPred (creating the python environment needed), you need to specify Scipion in the scipion.conf file both:

3.a) The path to the BepiPred folder as: BEPIPRED_HOME = <PathToBepiPred3_src> and

3.b) The activation command as: BEPIPRED_ACTIVATION_CMD = <ActivationCommand>

This way, Scipion will use your own BepiPred installation.


===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-chem-bepipred
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-bepipred.git

2. **Switch to the desired branch** (main or devel):

Scipion-chem-bepipred is constantly under development.
If you want a relatively older an more stable version, use main branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-bepipred
            git checkout devel

3. **Install**:

.. code-block::

            scipion installp -p path_to_scipion-chem-bepipred --devel

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg
