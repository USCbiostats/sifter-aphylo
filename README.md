
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

-----

The file

gene_association.goa_uniprot.80.txt which is formated in GAF 1.0
(see http://www-legacy.geneontology.org/GO.format.gaf-1_0.shtml)

Was downloaded from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/ version 80

The file hundredfamilies was downloaded from supplemental materials of the 2011
paper published in Genome Research. 

https://genome.cshlp.org/content/21/11/1969/suppl/DC1

https://ndownloader.figstatic.com/files/6649398 (see also https://doi.org/10.1371/journal.pcbi.0010045.sd001)
From the paper website https://doi.org/10.1371/journal.pcbi.0010045 (Original 2015 paper from plos comp bio)


 
The nudix family tree was downloaded from PFAM version 20 (as mentioned in the paper)
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam20.0/, and in particular, the
file trees.tar.gz


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

## Running SIFTER

To run sifter, in our case, we had to include the password to the databse. In particular,
following the example to "run SIFTER on some families":

```
python sifter_prepare.py --dbpass=[pass] -f PF12491,PF13820 path/to/families_data ../example/queries
python sifter_run.py ../example/queries ../example/results
python sifter_extract.py --dbpass=[pass] -f PF12491,PF13820 ../example/results ../examples/preds.txt
```

This also included creating the `../examples/` path which does not exists by default.

The **evidence files** are located in

`SIFTER-master/large_scale_v1.0/data/families_data/annotations/*.pli.gz`

## Manual

Under the hood, this is what is run

```
java -jar core/sifter2.1.1.jar
ERROR: Please specify run name:
usage: java -jar sifter.jar [OPTIONS] FAMILYNAME
 -sfx,--scale <filename>                     Set family .fx scale filename
                                             (default: data/scale-<FAMILY>.fx)
 -afx,--alpha <filename>                     Set family .fx alpha filename
                                             (default: data/alpha-<FAMILY>.fx)
 -fx,--familyfile <filename>                 Set family .fx parameter
                                             filename (default: data/infer-<FAMILY>.fx)
 -exp,--with-exp                             (Experiment) Use GOA protein
                                             annotations inferred from experiment.
 -floats,--with-floats                       Use floating numbers as
                                             confidence values.
 -genemania,--with-genemania                 Use protein predictions from
                                             genemania.
 -gene3d,--with-genethreed                   Use gene3d predictions in
                                             InterPro.
 -hamap,--with-hamap                         Use hamap predictions in
                                             InterPro.
 -iba,--with-iba                             (Computational) Use GOA
                                             protein annotations from those inferred from biological aspect of
                                             ancestor.
 -ibd,--with-ibd                             (Computational) Use GOA
                                             protein annotations from those inferred from biological aspect of
                                             descendant.
 -ic,--with-ic                               (Curator Statement) Use GOA
                                             protein annotations from those inferred from curator
 -ida,--with-ida                             (Experiment) Use GOA protein
                                             annotations inferred from direct assay.
 -iea,--with-iea                             (Automatically Assigned) Use
                                             GOA protein annotations inferred by electronic annotation.
 -iep,--with-iep                             (Experiment) Use GOA protein
                                             annotations from those inferred from expression profiles.
 -igc,--with-igc                             (Computational) Use GOA
                                             protein annotations from those inferred from genomic context.
 -igi,--with-igi                             (Experiment) Use GOA protein
                                             annotations from those inferred from genetic interaction.
 -ikr,--with-ikr                             (Computational) Use GOA
                                             protein annotations from those inferred from key residues.
 -imp,--with-imp                             (Experiment) Use GOA protein
                                             annotations from those inferred from mutant phenotype.
 -ipi,--with-ipi                             (Experiment) Use GOA protein
                                             annotations from those inferred from physical interaction.
 -ird,--with-ird                             (Computational) Use GOA
                                             protein annotations from those inferred from rapid divergence.
 -isa,--with-isa                             (Computational) Use GOA
                                             protein annotations from those inferred from sequence alignment.
 -ism,--with-ism                             (Computational) Use GOA
                                             protein annotations from those inferred from sequence model.
 -iso,--with-iso                             (Computational) Use GOA
                                             protein annotations from those inferred from sequence orthology.
 -iss,--with-iss                             (Computational) Use GOA
                                             protein annotations from those inferred from sequence similarity.
 -nd,--with-nd                               (Curator Statement) Use GOA
                                             protein annotations from those inferred from curator with no biological
                                             data available
 -nr,--with-nr                               (Obsolete) Use protein
                                             annotations with source not recorded.
 -folds,--folds <number>                     Number of folds in cross
                                             validation, leave-one-out is 0
 -pagosub,--with-pagosub                     Use protein predictions from
                                             pagosub.
 -panther,--with-panther                     Use panther predictions in
                                             InterPro.
 -pfam,--with-pfam                           Use Pfam predictions in
                                             InterPro.
 -pirsf,--with-pirsf                         Use PIRSF predictions in
                                             InterPro.
 -prints,--with-prints                       Use prints predictions in
                                             InterPro.
 -prodom,--with-prodom                       Use ProDom predictions in
                                             InterPro.
 -prosite_patterns,--with-prosite-profiles   Use PROSITE profile
                                             predictionss in InterPro.
 -protfun,--with-protfun                     Use protein predictions from
                                             ProtFun.
 -rca,--with-rca                             (Computational) Use GOA
                                             protein annotations from those reconstructed from computational analyses.
 -smart,--with-smart                         Use SMART predicitons in
                                             InterPro.
 -step,--step <number>                       Step size for gradient ascent
                                             in EM (M-step)
 -cutoff,--cutoff <number>                   Cutoff delta for gradient
                                             ascent in EM (M-step)
 -superfamily,--with-superfamily             Use SUPERFAMILY predictions
                                             in InterPro.
 -tas,--with-tas                             (Author Statement) Use GOA
                                             protein annotations from traceable author statements.
 -tigrfams,--with-tigr                       Use TIGRFAMs predictions in
                                             InterPro.
 -nas,--with-nas                             (Author Statement) Use GOA
                                             protein annotations from non-traceable author statements.
 -xval,--xvalidation                         (Run mode) Use
                                             cross-validation with EM.
    --help                                   Show help for arguments.
                                             (More help is available via README.txt)
 -em,--em                                    (Run mode) Perform EM to
                                             estimate parameters
 -g,--generate                               (Run mode) Generates a set of
                                             input parameters for the inference problem.
 -iter,--iter <number>                       Number of iterations. At the
                                             moment, this applies only to EM.
 -ontology,--ontology <filename>             Specify which ontology file
                                             you want (default: "data/function.ontology")
 -output,--output <filename>                 Set output file (default:
                                             output_directory/default.rdata)
 -phylo,--reconciled <filename>              Set reconciled .xml tree
                                             (default: reconciled/reconciled_<FAMILY>.xml)
 -pli,--protein <filename>                   Set protein file (default:
                                             proteins/proteinfamily_<FAMILY>.pli)
 -truncation,--truncation <number>           Number of functions to
                                             truncate to in approximation
 -v,--verbose                                Verbose operation.
```
