



# Instuctions for running Firedrake and access to "Anna's" desktop with multiple cores

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
Notes for old 2011 iMac with OS High Sierra 10.13.6
Login as toor:
```
>>> su toor
```
then type (give su-password):
```
>>> sudo python3.8 -m pip install numpy
```


