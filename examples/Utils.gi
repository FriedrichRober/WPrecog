
BindGlobal("RandomGenerator_ExtendGens", function(gens)
    local n, A, B;
    n := Length(gens);
    A := List([1 .. n], i -> Random([-1,1]));
    B := QuickRandomPermutedList(n);
    Add(gens, Product([1 .. n], i -> gens[B[i]]^A[i]));
end);

BindGlobal("RandomGenerator_SingleReplacement", function(gens, args...)
    local n, A, B, F, g, p;
    n := Length(gens);
    A := List([1 .. n], i -> Random([-1,0,1]));
    if Length(args) = 0 then
        F := Filtered([1 .. n], i -> A[i] <> 0);
        if Length(F) = 0 then
            return;
        fi;
        p := Random(F);
    else
        p := args[1];
        A[p] := 1;
    fi;
    B := QuickRandomPermutedList(n);
    gens[B[p]] := Product([1 .. n], i -> gens[B[i]]^A[B[i]]);
end);

BindGlobal("RandomGenerators", function(gens)
    local i, k;
    k := Length(gens);
    for i in [1 .. 5] do
        RandomGenerator_ExtendGens(gens);
    od;
    for i in [1 .. k] do
        RandomGenerator_SingleReplacement(gens, i);
    od;
    for i in [1 .. 10] do
        RandomGenerator_SingleReplacement(gens);
    od;
end);

BindGlobal("RandomGroup", function(G)
    local gens;
    gens := List(GeneratorsOfGroup(G));
    RandomGenerators(gens);
    return Group(gens);
end);
