LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Imprimitive
n := 9;
m := 12;
K := AlternatingGroup(n);
H := SymmetricGroup(m);

# Construct wreath product
PP := SymmetricGroup(n*m);
P := WreathProduct(SymmetricGroup(n), SymmetricGroup(m));
W := Group(Concatenation(
    List(GeneratorsOfGroup(K), x -> x^Embedding(P,1)),
    List(GeneratorsOfGroup(H), x -> x^Embedding(P,m+1))
));

# Random conjugation
c := PseudoRandom(PP);
G := W^(c^(-1));

# # Nice isomorphism
# isoWreath := IsomorphismWreathProduct(W);;
# isoConj := GroupHomomorphismByFunction(G, W, g -> g^c);
# isoNice := isoConj * isoWreath;;

# Recognition
data := RecogniseWreathProduct(RecogNode(G), rec(
    simpleBaseComp := false,
));
