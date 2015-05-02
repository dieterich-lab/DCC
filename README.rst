*****************************************
DCC: detect circRNAs from chimeric reads
*****************************************
DCC is a python package intended to detect and quantify circRNAs with high specificity. DCC works with the STAR chimeric.out.junction 
files which contains chimerically aligned reads including circRNA junction spanning reads. DCC detect and quantify circRNA junction 
spanning reads with high reliability. 

=================
Installation
=================
1) From git clone

.. code-block:: bash

  $ git clone git@github.com:dieterich-lab/DCC.git
  $ cd DCC
  $ python setup.py install
  
