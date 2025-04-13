LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Matrix Group : Imprimitive
d := 9;
q := 5;
K := GL(d,q);
m := 12;
H := Random(AllTransitiveGroups(NrMovedPoints, m));
D := d*m;

# Construct wreath product
P := GL(D, q);
W := WreathProduct(K, H);

# Random conjugation
c := PseudoRandom(P);;
G := W^(c^(-1));;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    assumeSimpleBaseComponent := false,
    recogniseBaseComponentBeforeDomain := false,
));;
output.res;

B := output.data.B;
Display(G.1^B);