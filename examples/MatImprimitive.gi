Read("gap/Experiment.gi");

EmbedGenerator := function(mat, m, q)
    local n, big, i, j;
    n := Length(mat);
    big := IdentityMat(m, GF(q));
    for i in [1..n] do
        for j in [1..n] do
            big[i][j] := mat[i][j];
        od;
    od;
    return big;
end;

# Matrix Groups : Imprimitive
d := 9;
q := 5;
K := GL(d,q);
m := 12;
H := SymmetricGroup(m);
# H := Random(AllTransitiveGroups(NrMovedPoints, m));
W := WreathProduct(K, H);
# isoWreath := IsomorphismWreathProduct(W);
D := d*m + 2;
P := GL(D, q);
W2 := Group(List(GeneratorsOfGroup(W), g -> EmbedGenerator(g, D, q)));
c := PseudoRandom(P);
G := W2^(c^(-1));
# conj := GroupHomomorphismByFunction(G, W, g -> g^c);
# iso := conj * isoWreath;;

BlocksMat := function(G)
    S := BlindDescend(G);
    # IsSingleComponent(S, iso);

    # S acts on W + U
    # W invariant space, dim(W) = d
    # U fixed space, dim(U) = (m-1)*d
    M := GModuleByMats(GeneratorsOfGroup(S), GF(q));
    v1 := RandomVector(D,q);
    v2 := RandomVector(D,q);
    # basis of W + < [v]\rho_U >
    L1 := MTX.SubmoduleGModule(M, v1);;
    L2 := MTX.SubmoduleGModule(M, v2);;
    # Length(L1) = d + 1;
    # Length(L2) = d + 1;
    W := Intersection(VectorSpace(GF(q), L1), VectorSpace(GF(q), L2));
    # Dimension(W) = d;
    tmp := TopGroupDomainVectorSpace(G, W);
    return tmp;
end;

tmp := BlocksMat(G); time;
domain := tmp.domain;
transversal := tmp.transversal;
B := Concatenation(List(domain, Basis));
# B1 := Concatenation(List(domain, Basis));
# B2 := SemiEchelonBasis(VectorSpace(GF(q), B1));
# B3 := Subspace(GF(q)^D, List(Basis(GF(q)^D), v->SiftedVector(B2,v)));
# B := Concatenation(B1, Basis(B3));
Display(G.3^(B^(-1)));