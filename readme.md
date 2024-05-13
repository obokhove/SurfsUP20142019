Install Jupiter (20-11-2023 from EP): 

>> mkdir jupiter

>> python3 -m venv .venv

>> ls -al

>> ls -al .venv

>> source .venv/bin/activate

>> which python 

>> deactivate

>> which python3

>> source .venv/bin/activate

>> pip install jupyterlab

>> pip install pandas

>> jupyter lab

(Opens then in browser.)

## Installation reinstall Firedrake on new macbook 02-05-2024
updated .zshrc; udated Xcode; update python in brew; still issue is that pkgutil is not found.
Frustrating: now trying to enforece brew's Python to be used instead of the mac one (all very unclear):
/opt/homebrew/bin/python3 firedrake-install

Hack via Jose pending repair of firedrake-install on chaco:
comment out pkgs.add("chaco") from line 722 of the install script. Do not comment out the return pkgs statement.
Use the --with-parmetis argument when you run the install script.
After installation, change line 1259 of firedrake/mesh.py from partitioner_type = "chaco" to partitioner_type = "parmetis"



# Instructions for running Firedrake and access to "Anna's" desktop with its 56 multiple cores and my imac:

## Wave-energy device pages


## Installation on macbook 13-06-2023 to 22-06-2023

First install xcode when asked 
login as:

>> su toor

then type in command from homebrew site without sudo
then as su:

>> brew install python3

>> pip3 install matplotlib

>> pip3 install lmfit

>> pip3 install pandas

>> pip3 install --upgrade matplotlib

>> brew doctor

>> brew cleanup

>> brew install php

Logged in as normal user then:

>> curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install

>> python3 firedrake-install

fails No such file or directory : 'brew'

Solved via assists using chatgpt: change path in .zshrc (as user) and use sudo chown (etc) (in admin) to give brew rights to user.
Detailed commands later:

As user, nano .zschrc in home directory of user; add (for me) location of homnebrew:

>> export PATH="/opt/homebrew/bin:$PATH"

then give access righst to user (as admin, "amtob" is my user account name):

>> sudo chown -R amtob /opt/homebrew

Works it seems:

>> . /Users/amtob/firedrake/bin/activate
or
>> . /Users/onnobokhove/amtob/werk/firedrake/bin/activate


## GitHub stuff

16-04-2023: Paraview crashed. Restart, then: cd /Users/onnobokhove/amtob/werk/firedrake/src/firedrake
and source firedrake/bin/activate
in

copied Eagre locally;
git clone https://github.com/obokhove/EagreEUEID20202023.git

and then added new code;
failed to commit it back onto eagre github
failed incomprehensible BS instructions everywhere; incomplete junk.


## 11-11-2022
Trying to install Firedrake using commands for FD-page:
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
then
python3 firedrake-install
Failed on petsc package??? but error message unclear.

Then did:
brew install python3
3x; notices python3 -V yields Python 3.10.8 but then my own non-FD Python code failed packages so typed:
>>> python3 -m pip install matplotlib
>>> python3 -m pip install --upgrade matplotlib
>>> python3 -m pip install lmfit
>>> python3 -m pip install pandas
Then my (tax-code) example worked again.
Removed firedrake with rm -rf firedrake
Then again
python3 firedrake-install

Also again added:
python3 -m pip install matplotlib 
python3 -m pip install pandas 
python3 -m pip install lmfit

16-11 Multiple delete firedrake, and re-installs, now by completing commenting out .zprofile and .zshrc files per instruction to move out the paths.
Still faild at pytest but install succeeded. Now touched .zchrc and .zprofile (getting a bit tiring).

26-11 Contunued. Deleted and reinstalled. Firedrake seems installed but Junho's periodic program does not run. Do not understand how to run:
cd $VIRTUAL_ENV/src/firedrake
pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
Ran as:
cd /Users/onnobokhove/amtob/werk/firedrake/src/firedrake  
then
./pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
or
onnobokhove@Os-MacBook-Pro firedrake % pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"  
zsh: command not found: pytest

Tried and failed:
onnobokhove@Os-MacBook-Pro werk % . /Users/onnobokhove/amtob/werk/firedrake/bin/activate
(firedrake) onnobokhove@Os-MacBook-Pro werk % cd $VIRTUAL_ENV/src/firedrake
pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
ImportError while loading conftest '/Users/onnobokhove/amtob/werk/firedrake/src/firedrake/tests/conftest.py'.
tests/conftest.py:9: in <module>
    from firedrake.utils import complex_mode
firedrake/__init__.py:80: in <module>
    from firedrake.output import *
firedrake/output.py:13: in <module>
    from .paraview_reordering import vtk_lagrange_tet_reorder,\
firedrake/paraview_reordering.py:5: in <module>
    from vtkmodules.vtkCommonDataModel import (
E   ModuleNotFoundError: No module named 'vtkmodules.vtkCommonDataModel'

GitHub-Firedrake; told to do: pip uninstall vtk and then pip install vtk.
    
    python3 -V
    
Then ran:
    
python3 -m pip uninstall vtk
    
but was not installed; then ran:
    
python3 pip install vtk

and was advised to run
     python3.10 -m pip install --upgrade pip
    so did that. Then start again installing FD (from werk -directory):
rm -rf firedrake
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
python3 firedrake-install
-waiting.
    Still failed the test; now did python3.10 -m pip uninstall vtk and then python3.10 pip install vtk uninstalled firedrake; started again.
Have to activate Firedrake and uninstall and install vtk within activated Firedrake, so:
. /Users/onnobokhove/amtob/werk/firedrake/src/firedrake/firedrake/bin/activate
then:
    
python3.10 -m pip uninstall vtk 
python3.10 -m pip install vtk  

Passed test on Friday 25-11-2022; FD-programs now seem to work.

15-12-2022: start up FD using:
. /Users/onnobokhove/amtob/werk/firedrake/src/firedrake/firedrake/bin/activate
then type:
cd $VIRTUAL_ENV/src/firedrake
pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
 


## Remote access to Anna's desktop

Replace amtob by your username:
```
>>> ssh -Y amtob@remote-access.leeds.ac.uk
>>> ssh -X amtob@mat-lin4123
>>> firedrake
>>> python
```
or 
```
>>> ssh -Y -J amtob@remote-access.leeds.ac.uk amtob@mat-lin4123
```
Then run Firedrake after installation:
```
>>> source firedrake/bin/activate
>>> firedrake/bin/activate [means a directory where the activate file of your firedrake is] 
>>> python main.py 
```
Currently using cumbersome two-step sftp to access mat-lin4123:
```
>>> sftp amtob@remote-access.leeds.ac.uk
```
and put files in appropriate direction; then login 
```
>>> ssh -Y -J amtob@remote-access.leeds.ac.uk
```
and sftp from there (from right dirctory) as follows:
```
>>> sftp amtob@mat-lin4123
```
or
```
>>> scp host_b:/scratch/tmp1/obokhove/data/BLf/* /Users/bokhoveo/dropbox/Variationalwaterwavemodels/NumericaltankJune2018/BL
```
or
```
>>> scp host_b:/scratch/tmp1/obokhove/data/BLf/* /any/path/athome
```
with in config file in .ssh directory:
 Host host_a
  User amtob
  Hostname remote-access.leeds.ac.uk

Host host_b
  User amtob
  Hostname mat-lin4123
  Port 22
  ProxyCommand ssh -q -W %h:%p host_a

## Firedrake

Following does work (newest version is firedrake-new):
```
>>> cp ~/vuurdraak/junho/* /scratch/tmp1/obokhove
>>> firedrake
>>> python KP_sol.py
```
MPI (see also Firedrake page; add "from firedrake.petsc import PETSc" in top code and adjust "print"--commands); works for up to 16 cores, fails for 20 cores, likely best for 14 cores since 14 cores per socket?:
```
>>>mpiexec -n 10 python 3D_tank.py
```

## Remote access to imac desktop
Same for imac at office (note that it is in sleep-mode so it needs some time to wake up before a login is possible, i.e. may need to try a few times while waiting for imac to wake up) (does not work 08-04-2020; works after a reboot 17-06-2020):
```
>>> ssh -Y -J amtob@remote-access.leeds.ac.uk amtob@mat-mac4175
```
22-10-2021 replaced by, in one window:
ssh -L 8888:mat-mac4175.leeds.ac.uk:22 amtob@remote-access.leeds.ac.uk
then in another window:
ssh -p 8888 amtob@localhost




## Updating old 2011 imac
Notes for old 2011 iMac with OS High Sierra 10.13.6

Fix 28-04: terminal gave command not known remar; fixed per terminal by running:
```
>>> PATH="/usr/bin:/bin:/usr/sbin:/sbin"
```
but that must then be done per terminal; fixed then by adding export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin" to .profile file; therefater when opening a terminal matters seemed restored.

Login as toor:
```
>>> su toor
```
then type (give su-password after each sudo; subsequently test program worked):
```
>>> sudo python3.8 -m pip install --upgrade certificates
>>> sudo python3.8 -m pip install numpy
>>> sudo python3.8 -m pip install scipy
>>> sudo python3.8 -m pip install matplotlib
>>> sudo -H python3.8 -m pip install --upgrade matplotlib
>>> sudo python3.8 -m pip install lmfit
>>> sudo python3.8 -m pip install pandas
```
31-03-2021 from user terminal and also within Firedrake:
```
>>> pip3 install numpy 
```
and for the other packages as in the above.


## Old notes Anna's machine
Dear all,

A summary of what was discussed earlier about remote computing.

[1] To connect to my machine remotely, you can follow either of the following two routes, but I suggest the second one. 

I. You have to first connect to amsta and then to mat-lin4123 (the name of the machine I am currently using):

ssh -X matak@amsta.leeds.ac.uk

ssh -Y matak@mat-li4123

II. Alternatively, you can use local port forwarding and the command below, which forwards local port 7777 to remote port 22, through amsta:

ssh -L 7777:mat-lin4123.leeds.ac.uk:22 -X matak@amsta.leeds.ac.uk

After this, you THEN OPEN A NEW TAB IN THE TERMINAL. In the new tab, you can then ssh directly on mat-lin4123 by typing:

ssh -p 7777 -Y matak@127.0.0.1

[2] To securely copy files from your local PC to the remote machine (mat-lin4123):

scp -P 7777 /local/path/to/filename matak@127.0.0.1:/path/to/filename

To securely copy files from the remote machine (mat-lin4123) to your local PC:

scp -P 7777 matak@127.0.0.1:/path/to/filename /local/path/to/filename

These two commands need to be executed while you are connected to your local PC.

[3] To run firedrake or anything else, you write (for firedrake you have to first use the source from the right directory):

nohup python filename.py >& output.out &

Nohup comes from “no hung up” after logout and the & symbol is telling it to run in the background.
After this you will get the process ID (PID) of the job running.

Alternatively, you can run the attached script file:

./script

which will execute the above command in the background and also sent you an email when it’s finished.

[4] To visualise or post-process the results, type:

paraview results.pvd &
matlab33 &
vlc movie.avi &

etc.

[5] To see what’s running and their PIDs, type:
ps -u

To kill a process: 
kill -9 PID1 PID2 ...

I can confirm that all of these commands work for me. Let me know if something doesn’t work for you and we should be able to figure it out!


## SSH key on github
ls -al ~/.ssh
ssh-keygen -t rsa -b 4096 -C "amtob@leeds.ac.uk"

vi id_rsa.pub
copy into GitHub; see: https://help.github.com/en/enterprise/2.15/user/articles/adding-a-new-ssh-key-to-your-github-account

When it works, e.g use/copy location from top-right in directory one wishes to clone in GitHub, say:
git clone git@github.com:colinjcotter/obcc_waves.git

## To webpage places
Tunnel in one window:
```
>>> ssh -L 8888:amsta.leeds.ac.uk:22 amtob@remote-access.leeds.ac.uk
```

and then set up sftp n another window:
```
>>> sftp -P 8888 obokhove@localhost
```
22-10-2021:
You can setup the tunnel using:
ssh -L 8889:leo.leeds.ac.uk:22 amtob@remote-access.leeds.ac.uk
and then you would connect to leo using your usual password (same pwd!):
sftp -P 8889 amtob@localhost


On Homebrew:
When asked:  sudo chown -R $(whoami)
replace by  sudo chown -R localusername e.g.  sudo chown -R onnobokhove (per instructions after brew doctor).
Also use echo 'export PATH="/usr/local/opt/python@3.8/bin:$PATH"' >> ~/.zshrc ???







