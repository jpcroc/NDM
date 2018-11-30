This is CEA's version of Iradina with a Graphical User Interface (GUI) designed to do SRIM like calculations.
This development has been made by Jean-Paul Crocombette and Christian Van Wambeke (CEA Saclay, France), under GPL license https://www.gnu.org/licenses/gpl.html
Bug reports or proposition of improvements should be addressed to Jean-Paul Crocombette jpcrocombette@cea.fr


installation notes for iradiana_GUI on linux (file iradinaGUI*.tgz)

Copy iradinaGUI*.tgz on your disk. Uncompress the tgz.

The root of the doumentation is in 
$HOMEIRA/iradinaGUI/doc/build/html/index.html
$HOMEIRA is the directory where you have put the tgz. $HOMEIRA contains the "iradinaGUI" directory.In $HOMEIRA you must create an iradina workir directory, (e.g. IRADINAGUI_WORKDIR) in which a LOG directory must also be created 

As explained in $HOMEIRA/iradinaGUI/doc/build/html/installation.html
If miniconda is NOT installed on your system you must go through steps 1 to 3. Otherwise go directly to step 4
1/ download miniconda from https://conda.io/miniconda.html. This gives you a file named Miniconda3-latest-Linux-x86_64.sh

2/>bash Miniconda3-latest-Linux-x86_64.sh . This actually installs miniconda

3/add the minionda path to your path >export PATH=where_you_put_it/miniconda3/bin:$PATH

Miniconda is installed on your system.

4/>conda create --name py3qt5 python=3 \
   pip jupyter matplotlib numpy pandas pandas-datareader \
   pyqt=5 scipy sympy jsonschema pyyaml libxml2 paramiko
This builds "py3qt5" which will be the environment on which the GUI will be run

5/ The following variables should be declared either in your .bashrc or in a script
export HOMEIRA="where you put the GUI directory"
export IRADINA_ROOT_DIR=$HOMEIRA/iradinaGUI
export IRADINAGUI_WORKDIR=$HOMEIRA/IRADINAGUI_WORKDIR
export IRADINAGUI_LOGDIR=$HOMEIRA/IRADINAGUI_WORKDIR/LOGS
export CONDADIR=~/miniconda3/ #where you put miniconda3

*/To Launch iradinaGUI : from the $HOMEIRA/iradinaGUI directory type
>source activate py3qt5
>iradinaGUI -g
>source deactivate
*/ An example of launching script is (assuming bash) :

export HOMEIRA=/home/croc/BCA/GUI_Iradina/linux
export IRADINA_ROOT_DIR=$HOMEIRA/iradinaGUI
export IRADINAGUI_WORKDIR=$HOMEIRA/IRADINAGUI_WORKDIR
export IRADINAGUI_LOGDIR=$HOMEIRA/IRADINAGUI_WORKDIR/LOGS
export CONDADIR=~/miniconda3/
source activate py3qt5
cd iradinaGUI
iradinaGUI -g
cd ../
source deactivate

 
TESTS ARE AVAILABLE IN iradinaGUI/examples

In the GUI :
The blue button with a question mark launchs the html help for the GUI and CEA installation. Please read it.
The green button with the question mark opens the pdf help of the iradina code in the original version by C. Borschel and C. Ronning.
Contextual help pops up when the mouse is placed over a variable or a button.
Some menus are accesssible with a right click.

