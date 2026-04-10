
BindGlobal("SetupImageViaProductAction", function(ri, data, options)
    local riS, domainDataS, n, m;
    m := data.m;
    n := RootInt(NrMovedPoints(Grp(ri)), m);
    riS := RecogNode(Group(List(GeneratorsOfGroup(data.S), StripMemory)));
    domainDataS := BlackBoxOrbit(riS, 1, rec(
        useMemory := false,
        maximalBoundOnTopDegree := n));
    data.domainS := domainDataS.domain;
    data.transversalS := domainDataS.transversal;
    if n <> Length(data.transversalS) then
        return fail;
    fi;
    data.translationS := MappingPermListList(data.domainS, [1..n]);
    data.parentWreathProduct := WreathProduct(SymmetricGroup(n), SymmetricGroup(m));
    return true;
end);

BindGlobal("ImageViaProductAction", function(g, ri, data, options)
    local pi, m, domain, t, riS, domainDataS, orb, s, n, W, sigma, w, i, j, k, points, p, candidates1, x, L, l, candidates2, images;
    pi := data.projectionFunction(g);
	if pi = fail then
		return fail;
	fi;
    m := data.m;
    t := List(data.transversal, StripMemory);
    orb := data.domainS;
    s := data.transversalS;
    n := Length(s);
    W := data.parentWreathProduct;
    sigma := data.translationS;
    w := EmptyPlist(m + 1);
    # now t_i*g*t_j^(-1) acts on orb like z_i*g_i*z_j^(-1),
    # hence like the image under phi that we want to define.
    for i in [1 .. m] do
        j := i^pi;
        points := [];
        for k in [1 .. n] do
            p := s[k]^(t[i] * g * t[j]^(-1));
            points[k] := RestrictedPerm(p, orb)^sigma;
        od;
        # try to find x = 1^w_i
        # Frequency Analysis: check that x has n distinct image points
        candidates1 := [];
        for x in [1 .. n] do
            L := [];
            for k in [1 .. n] do
                l := x^points[k];
                if l in L then
                    break;
                else
                    Add(L, l);
                fi;
            od;
            if Length(L) = n then
                Add(candidates1, x);
            fi;
        od;
        # Cycle Analysis: check that x occurs in the same cycle types in s[k]^(t[i] * g * t[j]^(-1)) as in s[k]
        if Length(candidates1) = 0 then
            return fail;
        elif Length(candidates1) = 1 then
            x := candidates1[1];
        else
            candidates2 := [];
            for x in candidates2 do
                if ForAll([1..n], k -> CycleLength(points[k], [1..n], x) = CycleLength(s[k], [1..n], 1)) then
                    Add(candidates2, x);
                fi;
            od;
            if Length(candidates2) = 1 then
                x := candidates2[1];
            else
                return fail;
            fi;
        fi;
        # now we found x = 1^w_i. Thus x^points[k] = k^w_i.
        images := [];
        for k in [1 .. n] do
            images[k] := x^points[k];
        od;
        w[i] := PermList(images);
    od;
    Add(w, pi);
    return WreathProductElementList(W, w);
end);
