LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Product Action
m := 4;
d := 4;
q := 3;
K := PSL(d,q);
n := NrMovedPoints(K);
n^m;
H := Random(AllTransitiveGroups(NrMovedPoints, m));

# Construct wreath product
P := SymmetricGroup(n^m);
W := WreathProductProductAction(K, H);;

# Random conjugation
c := PermList(QuickRandomPermutedList(NrMovedPoints(P)));;
G := W^(c^(-1));;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    action := "product action"
));;
res := output.res;
emb := output.data.embedding;;
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1));

timer := Runtime();;
riG := RecogniseGroup(G);;
timer := Runtime() - timer;;
Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds");