
################################
##                            ##
##          SL(n, q)          ##
##                            ##
################################

# parameters of classical group
n := 4;
p := 3;
e := 2;
q := p ^ e;

# classical groups
P := GL(n,q); # parent
M := SL(n,q); # matrix group in parent

# simple group, factoring scalar matrices out
ZP := Centre(P);
ZM := Centre(M);
quo := NaturalHomomorphismByNormalSubgroup(P, ZP);
T := Image(quo, M);

# automorphism group
# A = <I, D, G, F> = I ><| (<D, G> x F)
# I = T  and  F = C_e
# D = C_d, d = (n, q - 1)
# if n > 2  =>  G = C_2  and  <D, G> = D_{2 d}
# if n = 2  =>  G = C_1  and  <D, G> = C_{d}
A := AutomorphismGroup(T);
nice := NiceMonomorphism(A);

# inner automorphisms
Print("-------------------------------------------", "\n");
Print(" I  : Inner Automorphisms", "\n");
Print("|I| = ", Order(T), "\n");
emb_I := GroupHomomorphismByFunction(T, A, x -> ConjugatorAutomorphism(T, x));;
I := Image(emb_I); # = T

# diagonal automorphism
Print("-------------------------------------------", "\n");
Print(" D  : Diagonal Automorphisms", "\n");
d := Gcd(n, q - 1);
Print("|D| = ", d, "\n");
psi := NaturalHomomorphismByNormalSubgroup(Image(quo), T);
x := First(GeneratorsOfGroup(Image(psi)), x -> Order(x) = d);
y := PreImagesRepresentative(psi, x);
diag := y ^ (Order(y) / d);
diag_aut := ConjugatorAutomorphism(T, diag);
D := Group(diag_aut); # = C_d = F_q^* / F_q^{*n} = PGL(n, q) / PSL(n, q)
# C := Group(Z(q));
# lambda := GroupHomomorphismByFunction(C, C, c -> c ^ n);
# Cn := Image(lambda);
# rho := NaturalHomomorphismByNormalSubgroup(C, Cn);

# graph automorphisms
Print("-------------------------------------------", "\n");
Print(" G  : Graph Automorphisms", "\n");
# trivial
if n = 2 then
    Print("|G| = ", 1, "\n");
    G := Group(One(A)); # = C_1
# duality automorphism
else
    Print("|G| = ", 2, "\n");
    duality_aut := GroupHomomorphismByFunction(T, T, x -> (TransposedMat(PreImagesRepresentative(quo, x) ^ -1)) ^ quo);
    G := Group(duality_aut); # = C_2
fi;

# field automorphisms
Print("-------------------------------------------", "\n");
Print(" F  : Field Automorphisms", "\n");
Print("|F| = ", e, "\n");
frobenius := FrobeniusAutomorphism(GF(q));
frobenius_aut := GroupHomomorphismByFunction(T, T, x -> (PreImagesRepresentative(quo, x) ^ frobenius) ^ quo);
F := Group(frobenius_aut); # = C_e

# random element generator in A
sample_size := 100;
data := rec(
    emb_I := emb_I,
    frobenius_aut := frobenius_aut,
    diag_aut := diag_aut,
    G := G,
);;
opts := rec(

);;
RandomA := function(data, opts)
    local emb_I, frobenius_aut, diag_aut, G, T,
          elm_I, elm_F, elm_D, elm_G, exp_F, exp_D, exp_G, gen_G;
    emb_I := data.emb_I;
    frobenius_aut := data.frobenius_aut;
    diag_aut := data.diag_aut;
    G := data.G;
    # I
    T := Source(emb_I);
    elm_I := PseudoRandom(T) ^ emb_I;
    # F
    if IsBound(opts.exp_F) then
        exp_F := opts.exp_F;
    else
        exp_F := PseudoRandom([0 .. e - 1]);
    fi;
    elm_F := frobenius_aut ^ exp_F;
    # D
    if IsBound(opts.exp_D) then
        exp_D := opts.exp_D;
    else
        exp_D := PseudoRandom([0 .. d - 1]);
    fi;
    elm_D := diag_aut ^ exp_D;
    # G
    gen_G := GeneratorsOfGroup(G)[1];
    if IsBound(opts.exp_G) then
        exp_G := opts.exp_G;
    elif n = 2 then
        exp_G := 0;
    else
        exp_G := PseudoRandom([0, 1]);
    fi;
    elm_G := gen_G ^ exp_G;

    return rec(
        elm := elm_I * elm_F * elm_D * elm_G,
        elm_I := elm_I,
        exp_F := exp_F,
        elm_F := elm_F,
        exp_D := exp_D,
        elm_D := elm_D,
        exp_G := exp_G,
        elm_G := elm_G,
    );
end;;

sample := EmptyPlist(sample_size);
for i in [1 .. sample_size] do
    r := RandomA(data, opts);
    r.order := Order(r.elm);
    Add(sample, r);
    Print(r.order, ", ");
od;