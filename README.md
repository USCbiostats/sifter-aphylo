
# Running sifter on the aphylo trees

This repository contains the code used to make predictions with SIFTER using
the same phylogenetic trees and GO annotations featured in the aphylo paper.

## Downloads

The sifter program was downloaded from the GitHub repository [here](https://github.com/BrennerLab/SIFTER/tree/8274cac492a1c7355c4d92df548d6997dbb6b4b1).
The entire repository, which can be directly downloaded [here](https://github.com/BrennerLab/SIFTER/archive/8274cac492a1c7355c4d92df548d6997dbb6b4b1.zip), was decompress under `SIFTER-master`.

The data to run SIFTER was downloaded from

- SQL database https://sifter.berkeley.edu/media/sifter_db.gz
- Pfam families https://sifter.berkeley.edu/media/families_data.tar.gz

The Phylogenetic trees and GO annotations used to generate the data
from the aphylo paper can be obtained from here:

- Trees
- GO annotations

The current version of Biopython package was downloaded from [here](http://biopython.org/DIST/biopython-1.76.tar.gz)
as the version 1.76 is the latest to support Python 2.7, which is the version
used in SIFTER (we believe).

We needed to install pysqlite instead of sqlite, which needed the following lib:

```
sudo yum install libsqlite3x-devel.x86_64
```

Although it seems that that library is actually not imported anywhere in the program.

## Installing the system

1. Install MySQL and make sure there's enough space in your system to load the db.
   To do this, you need to change the `my.cnf` file that sets up where the db are
   stored. This can result in multiple problems in CentOS, including `SELinux`,
   which by default only allows `mysqld` to modify `/var/lib/mysql`.

2. In CentOS, MySQL-python can be installed directly from `yum`. 

3. The problem with 1. is that it results in errors in the case of 2., so I needed
   to create a symbolic link in the form of:

   ```
   sudo cd /var/lib && sudo ln -s /path/to/the/data/mysql mysql
   ```

4. Among the steps, need to use this
   
   ```
   ALTER USER 'yourusername'@'localhost' IDENTIFIED WITH mysql_native_password BY 'youpassword';
   ```
   
   Otherwise, the system returns with an error (source: [here](https://stackoverflow.com/questions/49194719/authentication-plugin-caching-sha2-password-cannot-be-loaded)).


