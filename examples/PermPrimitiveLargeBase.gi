LoadPackage("WPE", false);
LoadPackage("WPrecog", false);
ReadPackage("WPrecog", "examples/Utils.gi");

# Permutation Group : Product Action
m := 5;
l := 15;
k := 1;
l > 2*m*k^2; # Is Large Base?
n := Binomial(l, k);
n^m;
K := Action(AlternatingGroup(l), Combinations([1..l],k), OnSets);
H := Random(AllTransitiveGroups(NrMovedPoints, m));

# Construct wreath product
P := SymmetricGroup(n^m);
W := WreathProductProductAction(K, H);;
W := RandomGroup(W);;

# Random conjugation
c := PseudoRandom(P);;
G := W^(c^(-1));;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    action := "product action",
));;
res := output.res;
emb := output.data.embedding;;
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1));

timer := Runtime();;
riG := RecogniseGroup(G);;
timer := Runtime() - timer;;
Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds");