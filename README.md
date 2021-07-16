REQUIRED PACKAGES: numpy, pickle, jellyfish 

HOW TO USE:

Use control.py as the main file to control the program. Refer to control.py help (control.py -h) to see available options.

Example control.py terminal calls:

Perform all the steps and save the results (virus.txt, host.txt, dirs with pickle files and manhattan distance file) to ./k6 directory:
``` console
foo@bar:~$ python3 control.py -k 6 -o ./k6 -b ./host_dir/ -v ./virus_dir -p -d manhattan 
```

Prepare the files needed for distance calculations:
``` console
foo@bar:~$ python3 control.py -k 7 -o ./k7 -b ./host_dir/ -v ./virus_dir -p
```

Calculate Canberra distance on prepared files (requires the files from the picklify step):
``` console
foo@bar:~$ python3 control.py -o ./k7 -d canberra
```
kmerize.py, distance.py and picklify.py can also be used on their own if needed. 
