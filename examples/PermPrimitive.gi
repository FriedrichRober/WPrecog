LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Product Action
m := 5;
l := 15;
k := 1;
n := k * l;
l > 2*m*k^2; # Is Large Base?
K := Action(AlternatingGroup(n), Combinations([1..n],k), OnSets);
H := Random(AllTransitiveGroups(NrMovedPoints, m));

# Construct wreath product
PP := SymmetricGroup(n^m);
P := WreathProductProductAction(SymmetricGroup(n), SymmetricGroup(m));
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
    action := "product action",
));

RecogniseGroup(G);