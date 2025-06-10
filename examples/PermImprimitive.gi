LoadPackage("WPE", false);
LoadPackage("WPrecog", false);
ReadPackage("WPrecog", "examples/Utils.gi");

# Permutation Group : Imprimitive
n := 15;
m := 5;
K := AlternatingGroup(n);
H := SymmetricGroup(m);

# Construct wreath product
P := SymmetricGroup(n*m);
W := WreathProduct(K, H);;
W := RandomGroup(W);;

# Random conjugation
c := PseudoRandom(P);;
G := W^(c^(-1));;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    assumeSimpleBaseComponent := true,
));;
res := output.res;
emb := output.data.embedding;;
R := Range(emb);;
ListWreathProductElement(R, ImageElm(emb, G.1));