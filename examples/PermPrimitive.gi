LoadPackage("WPE", false);
LoadPackage("WPrecog", false);
ReadPackage("WPrecog", "examples/Utils.gi");

# Permutation Group : Product Action
m := 4;
d := 4;
q := 3;
K := PSL(d,q);
n := NrMovedPoints(K);
n^m; # 2 560 000
H := Random(AllTransitiveGroups(NrMovedPoints, m));

# Construct wreath product
P := SymmetricGroup(n^m);
W := WreathProductProductAction(K, H);;
W := RandomGroup(W);;

# Random conjugation
c := PermList(QuickRandomPermutedList(NrMovedPoints(P)));;
G := W^(c^(-1));;

# Recognition
emb := WreathProductDecomposition(RecogNode(G), rec(
    action := "product action"
));; # time ~ 3 min 20 seconds
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1)); # time ~ 10 seconds

timer := Runtime();;
riG := RecogniseGroup(G);;
timer := Runtime() - timer;;
Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds"); # time ~ 