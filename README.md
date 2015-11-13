# lattice_data

In this repository we will store the raw data for the [Lattices database](https://github.com/LMFDB/lmfdb-inventory/blob/master/db-Lattices.md) for [LMFDB](https://github.com/LMFDB/lmfdb).

Each line in a data file correspond to one integral lattice and it is a list containing the following items: 

* dimension, (integer)
* determinant, (integer)
* level, (integer)
* Gram matrix, (list of rows with integer entries) 
* density, (string)
* hermite number, (string)
* minimum, (integer)
* kissing number, (integer)
* [list of shortest vectors],  (list of lists of integers)
* size of aut. group, (integer)
* [theta series (coefficients)], (list of integers)
* class number, (integer)
* [list of matrices for genus representatives], (list of lists rows with integer entries)
* "name", (string)
* "comments", (string)

If a field os unknown, it will be inserted as "". 
For example, this is a line corresponding to a lattice:

[2,114,228,[[5,1],[1,23]],"0.367796388140185506419151585313","0.468292905790846982930329998491",5,2,[[1,0]],2,"",1,[[[5,1],[1,23]]],"",""]

In this example, the fields corresponding to the theta series coefficients, name and comment are not known (meaning not computed by the person who upload the data) and so are inserted as "". 

The order of the entries is fixed: if the data is given in a different order, this needs to be recorded before inserting this in the LMFDB database.

This Repository contains the following files:

1. compiled-data: This contains the data on class number 1 lattices, it was uploaded on November 11th.

2. raw-data: This is the ascii file from the original Lorch and Kirschmer paper.

3. raw-data-list: This takes the ascii file and turns each line into a list item to be ready by Magma. 
