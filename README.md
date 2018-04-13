# Lattice data

In this repository we will store the raw data for the [Lattices database](https://github.com/LMFDB/lmfdb-inventory/blob/master/db-Lattices.md) for [LMFDB](https://github.com/LMFDB/lmfdb).

Each line in a data file correspond to one integral lattice and it is a list containing the following items: 

* dimension, (integer)
* negative_part_of_signature, (integer)
* determinant, (integer)
* level, (integer)
* Gram matrix, (list of rows with integer entries) 
* density, (string)
* hermite number, (string)
* minimum, (integer)
* kissing number, (integer)
* [list of shortest vectors],  (list of lists of integers)
* [size of aut. group, number of generator, [list of generators as matrices], group structure], (list of [int, int, [lists of ints], string, string])
* class number, (integer)
* [list of matrices for genus representatives], (list of lists rows with integer entries)
* "name", (string)
* "comments", (string)

If a field os unknown, it will be inserted as "". 
For example, this is a line corresponding to a lattice:
```
[2,114,228,[[5,1],[1,23]],"0.367796388140185506419151585313","0.468292905790846982930329998491",5,2,[[1,0]],2,"",1,[[[5,1],[1,23]]],"",""]
```
In this example, the fields corresponding to the theta series coefficients, name and comment are not known (meaning not computed by the person who upload the data) and so are inserted as "". 

The order of the entries is fixed: if the data is given in a different order, this needs to be recorded before inserting this in the LMFDB database.

The labeing scheme for lattices will be as follows,

***dimension.negative_part_of_signature.determinant.level.genus_enumerator.class_enumerator***

where 

* negative_part_of_signature will most likely be 0 since at present the database only contains positive definite integral lattices.
* genus_enumerator counts over the genera of a fixed dimension, determinant, and level. 
* class_enumerator counts over the classes of a fixed genus.




 
