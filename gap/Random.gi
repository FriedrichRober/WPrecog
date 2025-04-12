
# Returns a random permutation on n as a list
InstallGlobalFunction("QuickRandomPermList", function(n)
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

# Returns a random subproduct of G
InstallGlobalFunction("QuickRandomSubproduct", function(G)
    local gens, n, A, B;
    if IsGroup(G) then
        gens := GeneratorsOfGroup(G);
    else # G is a list of group elements
        gens := G;
    fi;
    n := Length(gens);
    A := List([1 .. n], i -> Random([-1,0,1]));
    B := QuickRandomPermList(n);
    return Product([1 .. n], i -> gens[B[i]]^A[i]);
end);

# normal closure of H in G, i.e. <H^G>
InstallGlobalFunction("QuickNormalClosure", function(H, G, n)
    local gens;
    if IsGroup(H) then
        gens := GeneratorsOfGroup(H);
    else
        gens := H;
    fi;
    return Group(FastNormalClosure(G, gens, n));
end);

BindGlobal("QuickRandomVector", function(d, q)
    local F;
    F := GF(q);
    return List([1..d], i -> PseudoRandom(F));
end);