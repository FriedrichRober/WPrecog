LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Product Action
m := 3;
d := 3;
q := 5;
K := PSL(d,q);
n := NrMovedPoints(K);
n^m;
H := Random(AllTransitiveGroups(NrMovedPoints, m));

# Construct wreath product
PP := SymmetricGroup(n^m);
P := WreathProductProductAction(K, H);;
W := P;;

# Random conjugation
c := PermList(QuickRandomPermutedList(NrMovedPoints(PP)));;
G := W^(c^(-1));;

# # Nice isomorphism
# isoWreath := IsomorphismWreathProduct(P);;
# isoConj := GroupHomomorphismByFunction(G, P, g -> g^c);
# isoNice := isoConj * isoWreath;;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    action := "product action",
));;
res := output.res;

RecogniseGroup(G);