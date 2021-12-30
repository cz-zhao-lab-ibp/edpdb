# EdPDB
EdPDB: A Multi-Functional Tool for Protein Structure Analysis.

__We recommend all the users run EdPDB in a container. Please refer to the [edpdb_docker](https://github.com/cz-zhao-lab-ibp/edpdb_docker) page for detailed information.__

EdPDB is a program to handle PDB (Brookhaven Protein Databank) format coordinate file(s). With EdPDB, a PDB file may work as a small database, where information like distances, angles etc. within the file can be easily obtained.

## Documentation
The documentation of EdPDB can be found in [this page](http://175.24.30.189:8081/edpdb/).

## Reference
Zhang, X-J., and Brian W. Matthews. "EDPDB: a multifunctional tool for protein structure analysis." Journal of applied crystallography 28.5 (1995): 624-630.

## The letter from author
Dear EDPDB users,

EDPDB program package consists of one program and a number of accessary files. An updated EDPDB homepage is kept at ~~http://www.ibp.cas.cn/zhangkLab/zhangklabedpdb/index.html.~~
This EDPDB distribution is a Linux "tar" file, "edpdb_linux_06a.tar". It was created using the following command.

```
"tar -cf edpdb_lx_v06a.tar edpdb"
```

where the directory "edpdb" contains some subdirectories which should be restored accordingly (using "tar -xf edpdb_lx_v06a.tar". To run the excutable "edpdb", one needs to include the file 
"edpdb.csh" in his/her ".login" file.  This "edpdb.csh" file defines the required logicals names, etc.

ONE MAY NEED TO MODIFY THIS FILE BEFORE USE IT ON YOUR OWN COMPUTER.

The excutable edpdb was compiled under Linux with the Fortran77 compiler-g77 (GNU project Fortran Compiler (v0.5.24)). 
The excutable usually can run on a Linux system. If the program has to be recompiled on your own computer, use the "makefile" in this package. You need to have the "GNU Readline" software package. In the "makefile", you should specify the location of the "GNU Readline" library using the flag id "LDFLAGS -L" and "CFLAGS -I".

I encourage you to report any bugs, and/or fixes if available.

With best regards,

X. Cai Zhang,  Ph.D.

## Contacts
© [X. Cai Zhang](http://english.ibp.cas.cn/faculty/index_18316.html?json=http://www.ibp.cas.cn/sourcedb_ibp_cas/cn/ibpexport/EN_xsszmZ/202005/t20200519_5582960.json). All Rights Reserved. 

Institute of Biophysics, CAS 

Email: zhangc@ibp.ac.cn 

Last modified on 2020-10-21 .

