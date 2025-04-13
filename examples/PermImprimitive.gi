LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Imprimitive
n := 9;
m := 12;
K := AlternatingGroup(n);
H := SymmetricGroup(m);

# Construct wreath product
P := SymmetricGroup(n*m);
W := WreathProduct(K, H);;

# Random conjugation
c := PseudoRandom(P);;
G := W^(c^(-1));;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    assumeSimpleBaseComponent := true,
));;
output.res;