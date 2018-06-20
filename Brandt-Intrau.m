/* 

This is magma code to compute LMFDB data for the Brandt-Intrau Latices available here (http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/Brandt_1.html).  The form

            ax^2+by^2+cz^2+dyz+exz+fxy
                        
has the associated matrix

            [2a f  e ]
        A = [f  2b d ]
            [e  d  2c]
                        
and determinant 

            8abc + 2fde - 2ad^2 - 2be^2 - 2cf^2 > 0
                        
since the form is positive definite.  Brandt and Intrau call these form odd if d,e,f are all even, and even otherwise. In the langauge of Conway and Sloane, the odd forms are type I_3 and even are of type II_3.  For odd forms, following the convetion for lattices in the LMFDB, we will give (1/2)A as the Gram matrix, and for even forms will will give A as the Gram matrix.

*/



SetColumns(0);

function ShortString(x)
s:=Sprint(x);
sx:="";
for i in [1..#s] do
if not s[i] in ["\t", "\n","*"] then
sx:=sx*s[i];
end if;
end for;
return sx;
end function;

/* The lists BIodd and BIeven contain lists of coefficients of odd and even forms, respectively. */

R:=MatrixRing(Integers(),3);

/* The code below appents the lists of coefficents in BIodd and BIeven to one central list which contains lists of coefficients for MAGMA lattices. */

LatticeCoeff:=[];
for v in BIeven do
Append(~LatticeCoeff,[2*v[1],v[6],2*v[2],v[5],v[4],2*v[3]]);
end for;

for v in BIodd do
Append(~LatticeCoeff,[v[1],(1/2)*v[6],v[2],(1/2)*v[5],(1/2)*v[4],v[3]]);
end for;

LatticeList:=[];
for v in LatticeCoeff do
Append(~LatticeList,LatticeWithGram(R!SymmetricMatrix(v)));
end for;

/* For testing purposes, I'm going to make a short list.  Comment this away when things are sorted.

ShortList:=[];
for i in [1..100] do
Append(~ShortList,LatticeList[i]);
end for;*/


/* The code below translates the lists of coefficents in BIodd and BIeven into coeffiencts for SAGE quadratic forms. 

FullList:=[* *];
for v in BIeven do
Append(~FullList,[2*v[1],2*v[6],2*v[5],2*v[2],2*v[4],2*v[3]]);
end for;

for v in BIodd do
Append(~FullList,[v[1],v[6],v[5],v[2],v[4],v[3]]);
end for;

Write("TernaryList",FullList);

*/

/* Now from LatticeList we can compute the data for the LMFDB as described here (https://github.com/annahaensch/lattice_data/blob/master/README.md). The data will write out to a file, line by line */

for V in LatticeList do
    DataList:=[* *];
    n:=Dimension(V);
    gram:=[];
    for i in [1..n] do
        row:=[];
        for j in [1..n] do
            Append(~row,GramMatrix(V)[i][j]);
        end for;
        Append(~gram,row);
    end for;

    S:=ShortestVectors(V);
    vec_list:=[];
    for s in S do
        vec:=[];
        for i in [1..n] do
            Append(~vec,s[i]);
        end for;
        Append(~vec_list,vec);
    end for;

    GenusRepresentatives:=[];
    for r in Representatives(Genus(V)) do
	if IsIsometric(r,V) eq true then 
	    mat:=[];
        for i in [1..n] do
            row:=[];
            for j in [1..n] do
                Append(~row,InnerProductMatrix(V)[i][j]);
            end for;
            Append(~mat,row);
        end for;
        Append(~GenusRepresentatives,ShortString(mat));
	else
        I:=Index(LatticeList,V);
        G:=#Genus(V);
        Low:=Max(1,I-G);
        High:=Min(#LatticeList,I+G);
        for w in LatticeList[Low..High] do
            if IsIsometric(r,w) eq true then
                mat:=[];
                for i in [1..n] do
                    row:=[];
                    for j in [1..n] do
                        Append(~row,InnerProductMatrix(w)[i][j]);
                    end for;
                    Append(~mat,row);
                end for;
                Append(~GenusRepresentatives,ShortString(mat));
                break;
            end if;
        end for;
    end if;
    end for;

    /* This is just a sanity check to see that we that we captured all of the genus representatives */
    if not #GenusRepresentatives eq #Genus(V) eq true then
        Append(~GenusRepresentatives,"incomplete");
    end if;
    /* It should always be complete for Brandt-Intrau since we know it is so.*/

    AutList:=[**];
    Append(~AutList,#AutomorphismGroup(V));
    Append(~AutList,#Generators(AutomorphismGroup(V)));
    GenList:=[];
    for v in Generators(AutomorphismGroup(V)) do
        GenMat:=[];
        for i in [1..n] do
            row:=[];
            for j in [1..n] do
                Append(~row,v[i][j]);
            end for;
            Append(~GenMat,row);
        end for;
        Append(~GenList,ShortString(GenMat));
    end for;
    Append(~AutList,ShortString(GenList));
    Append(~AutList,"\"\"");


    Append(~DataList,n);
    Append(~DataList,0);
    Append(~DataList,Determinant(V));
    Append(~DataList,Level(V));
    Append(~DataList,ShortString(gram));
    Append(~DataList,"\"" cat ShortString(Density(V)) cat "\"");
    Append(~DataList,"\"" cat ShortString(HermiteNumber(V)) cat "\"");
    Append(~DataList,Minimum(V));
    Append(~DataList,KissingNumber(V));
    Append(~DataList,vec_list);
    Append(~DataList,AutList);
    Append(~DataList,"\"\"");
    Append(~DataList,#Genus(V));
    Append(~DataList,ShortString(GenusRepresentatives));
    Append(~DataList,"\"\"");
    Append(~DataList,"This lattice appears in the Brandt-Intrau table of primitive positive definite ternary quadratic forms of discriminant up to 2000.");
    Write("data", ShortString(DataList));
    /* Write("data","");  uncomment this to get a line of blank space between entries.*/
end for;






