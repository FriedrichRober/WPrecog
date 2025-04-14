LoadPackage("recog", false);
LoadPackage("WPE", false);

DeclareGlobalFunction("BlindDescend");
DeclareGlobalFunction("QuickRandomList");
DeclareGlobalFunction("QuickRandom");
DeclareGlobalFunction("BlindDescend_Step");

InstallGlobalFunction("QuickRandomList", function(n)
    local lst, i, j, tmp;
    lst := [1..n];
    for i in [n, n-1..2] do
        j := Random([1..i]);
        tmp := lst[i];
        lst[i] := lst[j];
        lst[j] := tmp;
    od;
    return lst;
end);

InstallGlobalFunction("QuickRandom", function(gens)
    local n, A, B;
    n := Length(gens);
    A := List([1 .. n], i -> Random([-1,0,1]));
    B := QuickRandomList(n);
    return Product([1 .. n], i -> gens[B[i]]^A[i]);
end);

InstallGlobalFunction("BlindDescend_Step", function(gens, P, n)
    local y, ord, z;
    y := QuickRandom(gens);
    ord := Order(y);
    if IsEvenInt(ord) then
        z := y ^ (ord / 2);
        return FastNormalClosure(P, [z], n);
    else
        return gens;
    fi;
end);

InstallGlobalFunction("BlindDescend", function(args...)
    local G, l1, l2, n, n2, S, gens, i, B;

    if Length(args) < 1 then
        Error("Argument G is required.");
    fi;

    G := args[1];

    if Length(args) >= 2 then
        l1 := args[2];
    else
        l1 := 5;
    fi;

    if Length(args) >= 3 then
        l2 := args[3];
    else
        l2 := 35;
    fi;

    if Length(args) >= 4 then
        n1 := args[4];
    else
        n1 := 5; # roughly 30 generators
    fi;

    if Length(args) >= 5 then
        n2 := args[5];
    else
        n2 := 2;
    fi;

    S := G;
    gens := List(GeneratorsOfGroup(S));
    for i in [1 .. l1] do
        gens := BlindDescend_Step(gens, S, n1);
        S := Group(gens);
    od;
    B := S; # inside base group
    for i in [1 .. l2] do
        gens := BlindDescend_Step(gens, S, n1);
        S := Group(gens);
    od;
    # As a pre-caution, take normal closure in base group, if we descended too far
    S := Group(FastNormalClosure(B, gens, n2));
    return S;
end);

BindGlobal("IsSingleComponent", function(G, iso)
    local T, g, t;
    T := fail;
    for g in GeneratorsOfGroup(G) do
        t := Territory(g^iso);
        if Length(t) > 1 then
            return false;
        fi;
        if T = fail then
            T := t;
        elif T <> t then
            return false;
        fi;
    od;
    return true;
end);

# TODO, make faster for matrix group, by checking component-wise
DeclareGlobalFunction("DoCommuteElms");
InstallGlobalFunction(DoCommuteElms, function(x, y)
    return x * y = y * x;
end);

DeclareGlobalFunction("DoCommuteGroups");
InstallGlobalFunction(DoCommuteGroups, function(A, B)
    local a, b;
    for a in GeneratorsOfGroup(A) do
        for b in GeneratorsOfGroup(B) do
            if not DoCommuteElms(a, b) then
                return false;
            fi;
        od;
    od;
    return true;
end);

MergeGroups := function(Sig, Sj, Sk)
    return Group(Concatenation(
        GeneratorsOfGroup(Sk),
        GeneratorsOfGroup(Sig),
        GeneratorsOfGroup(Sj)
    ));
end;

MergeVectorSpaces := function(A, B, q)
    return VectorSpace(GF(q),
        Concatenation(
        Basis(A),
        Basis(B)
    ));
end;

AreVectorSpacesDisjoint := function(A, B)
    for a in Basis(A) do
        if a in B then
            return false;
        fi;
    od;
    return true;
end;

DeclareGlobalFunction("TopGroupDomain");
InstallGlobalFunction(TopGroupDomain, function(args...)
    local G, S, isDisjoint, merge, recursionDepth, recursionMaxDepth, maximalBoundOnTopDegree, t, domain, i, g, Si, Sig, Sj, foundNew, Sk, SS;

    if Length(args) < 2 then
        Error("Arguments G and S are required.");
    fi;

    G := args[1];
    S := args[2];

    if Length(args) >= 3 then
        recursionDepth := args[3];
    else
        recursionDepth := 1;
    fi;

    if Length(args) >= 4 then
        recursionMaxDepth := args[4];
    else
        recursionMaxDepth := 5;
    fi;

    if Length(args) >= 5 then
        maximalBoundOnTopDegree := args[5];
    else
        maximalBoundOnTopDegree := 50;
    fi;

    if recursionMaxDepth <= recursionDepth then
        return fail;
    fi;

    t := [One(G)];
    domain := [S];
    i := 1;

    while i <= Length(domain) do
        Si := domain[i];
        for g in GeneratorsOfGroup(G) do
            Sig := OnPoints(Si, g);
            foundNew := true;
            Sk := fail;
            for Sj in domain do
                if not DoCommuteGroups(Sig, Sj) then
                    foundNew := false;
                    if Sk = fail then
                        Sk := Sj;
                    else
                        SS := MergeGroups(Sig, Sj, Sk);
                        return TopGroupDomain(G, SS, recursionDepth + 1, recursionMaxDepth);
                    fi;
                fi;
            od;
            if foundNew then
                if Length(domain) >= maximalBoundOnTopDegree then
                    return fail;
                fi;
                Add(t, t[i] * g);
                Add(domain, Sig);
            fi;
        od;
        i := i + 1;
    od;

    return rec(domain := domain, transversal := t);
end);

DeclareGlobalFunction("TopGroupDomainVectorSpace");
InstallGlobalFunction(TopGroupDomainVectorSpace, function(args...)
    local G, S, q, recursionDepth, recursionMaxDepth, maximalBoundOnTopDegree, t, domain, i, g, Si, Sig, Sj, foundNew, Sk, SS;

    if Length(args) < 2 then
        Error("Arguments G and S are required.");
    fi;

    G := args[1];
    q := Size(DefaultFieldOfMatrixGroup(G));
    S := args[2];

    if Length(args) >= 3 then
        recursionDepth := args[3];
    else
        recursionDepth := 1;
    fi;

    if Length(args) >= 4 then
        recursionMaxDepth := args[4];
    else
        recursionMaxDepth := 5;
    fi;

    if Length(args) >= 5 then
        maximalBoundOnTopDegree := args[5];
    else
        maximalBoundOnTopDegree := 50;
    fi;

    if recursionMaxDepth <= recursionDepth then
        return fail;
    fi;

    t := [One(G)];
    domain := [S];
    i := 1;

    while i <= Length(domain) do
        Si := domain[i];
        for g in GeneratorsOfGroup(G) do
            Sig := OnPoints(Si, g);
            foundNew := true;
            for Sj in domain do
                if Sig = Sj then
                    foundNew := false;
                    break;
                elif not IsTrivial(Intersection(Sig, Sj)) then
                    SS := MergeVectorSpaces(Sig, Sj, q);
                    Basis(SS);
                    Dimension(SS);
                    return TopGroupDomainVectorSpace(G, SS, recursionDepth + 1, recursionMaxDepth);
                fi;
            od;
            if foundNew then
                if Length(domain) >= maximalBoundOnTopDegree then
                    return fail;
                fi;
                Add(t, t[i] * g);
                Add(domain, Sig);
            fi;
        od;
        i := i + 1;
    od;

    return rec(domain := domain, transversal := t);
end);

BindGlobal("TopComponentImage", function(g, domain)
    local m, images, i, j, Si, Sj, Sig;
    m := Length(domain);
    images := EmptyPlist(m);
    for i in [1 .. m] do
        Si := domain[i];
        Sig := OnPoints(Si, g);
        for j in [1 .. m] do
            if j in images then
                continue;
            fi;
            Sj := domain[j];
            if not DoCommuteGroups(Sig, Sj) then
                images[i] := j;
                break;
            fi;
        od;
        if not IsBound(images[i]) then
            return fail;
        fi;
    od;
    return PermList(images);
end);


BindGlobal("RandomVector", function(d, q)
    local F;
    F := GF(q);
    return List([1..d], i -> PseudoRandom(F));
end);

BindGlobal("StandardVector", function(d, q, i)
    local F, v;
    F := GF(q);
    v := ListWithIdenticalEntries(d, Zero(F));
    v[i] := One(F);
    return v;
end);

BindGlobal("NormaliseVector", function(d, q, b, I)
    local F, v, i;
    F := GF(q);
    v := ListWithIdenticalEntries(d, Zero(F));
    for i in I do
        v[i] := b[i];
    od;
    return v;
end);

BindGlobal("ExtendAndNormaliseBasis", function(d, q, L)
    local B, I, j, i;
    B := [];
    I := [];
    for j in [1 .. Length(L)] do
        Add(I, First([1 .. d], i -> L[j][i] = One(GF(q))));
    od;
    j := 1;
    for i in [1 .. d] do
        if j <= Length(L) and L[j][i] = One(GF(q)) then
            Add(B, L[j]);
            j := j + 1;
        else
            Add(B, StandardVector(d, q, i));
        fi;
    od;
    return B;
end);