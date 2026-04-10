LoadPackage("WPE", false);
LoadPackage("WPrecog", false);
ReadPackage("WPrecog", "examples/Utils.gi");

# Permutation Group : Product Action
m := 5;
l := 15;
k := 1;
l > 2*m*k^2; # Is Large Base?
n := Binomial(l, k);
n^m; # 759 375
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
emb := WreathProductDecomposition(RecogNode(G), rec(
    action := "product action",
));; # time ~ 10 seconds
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1)); # time ~ 1 seconds

# timer := Runtime();;
# riG := RecogniseGroup(G);;
# timer := Runtime() - timer;;
# Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds"); # time ~ 100 seconds

riG := FindHomMethodsPerm.LargeBasePrimitive(RecogNode(G), G);; # time ~ 100 seconds

# Permutation Group : Product Action
m := 7;
l := 8;
k := 1;
l > 2*m*k^2; # Is Large Base?
n := Binomial(l, k);
n^m; # 2 097 152
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
emb := WreathProductDecomposition(RecogNode(G), rec(
    action := "product action",
));; # time ~ 1 min
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1)); # time ~ 4 seconds

# timer := Runtime();;
# riG := RecogniseGroup(G);;
# timer := Runtime() - timer;;
# Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds"); # time ~ 

riG := FindHomMethodsPerm.LargeBasePrimitive(RecogNode(G), G);; # time ~ 11 min and fail