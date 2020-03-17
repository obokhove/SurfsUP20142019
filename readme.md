



# Instuctions for running Firedrake and access to "Anna's" desktop with multiple cores

## Remote access to Anna's desktop
```
>>> ssh -X username@mat-lin4123
>>> firedrake
>>> python
```

Make a tunnel (? not tested) in one terminal:
```
>>> ssh -L 8888:amsta.leeds.ac.uk:22 amtob@remote-access.leeds.ac.uk
```


## Run Firedrake after installation:

```
>>> source firedrake/bin/activate
>>> firedrake/bin/activate [means a directory where the activate file of your firedrake is] 
>>> python main.py 
```




