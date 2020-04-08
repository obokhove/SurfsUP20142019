



# Instructions for running Firedrake and access to "Anna's" desktop with multiple cores

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

Same for imac at office (note that it is in sleep-mode so it needs some time to wake up before a login is possible, i.e. may need to try a few times while waiting for imac to wake up):
```
>>> ssh -Y -J amtob@remote-access.leeds.ac.uk amtob@mat-mac4175
```


II. Notes for old 2011 iMac with OS High Sierra 10.13.6
Login as toor:
```
>>> su toor
```
then type (give su-password after each sudo; subsequently test program worked):
```
>>> sudo python3.8 -m pip install numpy
>>> sudo python3.8 -m pip install scipy
>>> sudo python3.8 -m pip install matlibplot
>>> sudo -H python3.8 -m pip install --upgrade matplotlib
>>> sudo python3.8 -m pip install lmfit
```


