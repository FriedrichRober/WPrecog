BindGlobal("DoCommuteElms", function(x, y)
    return x * y = y * x;
end);

BindGlobal("DoCommuteGroups", function(A, B)
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

BindGlobal("TopGroupDomainViaSimpleSubgroup", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: TopGroupDomainViaSimpleSubgroup(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> not DoCommuteGroups(a,b); # if groups do not commute, groups must be equal
    userOptions.isDistinguishable := {a, b} -> true; # if groups commute, groups must be disjoint

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("TopGroupDomainViaGenericSubgroup", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: TopGroupDomainViaGenericSubgroup(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> false; # we cannot check for equality without further assumptions
    userOptions.isDistinguishable := {a, b} -> DoCommuteGroups(a,b); # if groups commute, then they might lie in different components
    userOptions.isConflict := {p, P} -> Length(P) >= 2; # we found two groups with which p does not commute
    userOptions.resolveConflict := {p, P} -> Group(Concatenation(List(Concatenation(P, [p]), GeneratorsOfGroup))); # we merge these groups and hope that we now cover the socle in the component

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("TopGroupPointViaLinearAction", function(S)
    local d, q, M, v1, v2, L1, L2, W;
    d := DimensionOfMatrixGroup(S);
    q := Size(DefaultFieldOfMatrixGroup(S));
    # S acts on W + U
    # W invariant space, dim(W) = d
    # U fixed space, dim(U) = (m-1)*d
    M := GModuleByMats(GeneratorsOfGroup(S), GF(q));
    v1 := QuickRandomVector(d,q);
    v2 := QuickRandomVector(d,q);
    # basis of W + < [v]\rho_U >
    L1 := MTX.SubmoduleGModule(M, v1);;
    L2 := MTX.SubmoduleGModule(M, v2);;
    W := Intersection(VectorSpace(GF(q), L1), VectorSpace(GF(q), L2));
    return W;
end);

BindGlobal("TopGroupDomainViaLinearAction", function(args...)
    local ri, s, q, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: TopGroupDomainViaLinearAction(ri, s[, options])");
    fi;

    ri := args[1];
    q := Size(DefaultFieldOfMatrixGroup(Grp(ri)));
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> a = b;
    userOptions.isDistinguishable := {a, b} -> IsTrivial(Intersection(a, b));
    userOptions.isConflict := {p, P} -> Length(P) >= 1;
    userOptions.resolveConflict := function(p1, P)
        local p2, pNew;
        p2 := P[1];
        pNew := VectorSpace(GF(q), Concatenation(Basis(p1),Basis(p2)));
        Basis(pNew);
        Dimension(pNew);
        return pNew;
    end;

    return BlackBoxOrbit(ri, s, userOptions);

end);

BindGlobal("TopGroupDomainViaProductAction", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: TopGroupDomainViaProductAction(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> Set(Orbit(a, 1)) = Set(Orbit(b, 1));
    userOptions.isDistinguishable := {a, b} -> Size(Intersection(Orbit(a, 1), Orbit(b, 1))) = 1;
    userOptions.isConflict := {p, P} -> Length(P) >= 1;
    userOptions.resolveConflict := {p, P} -> Group(Concatenation(List(Concatenation(P, [p]), GeneratorsOfGroup)));

    return BlackBoxOrbit(ri, s, userOptions);
end);

InstallGlobalFunction("BlackBoxOrbit", function(args...)
    local ri, s, userOptions, options, name, G, t, domain, i, a, b, g, p, isEq, foundNew, indistinguishablePoints, j;

    # ==================================================
    # Input extraction
    # ==================================================
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: BlackBoxOrbit(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;

    # ==================================================
    # Default options
    # ==================================================
    options := rec(
        # -- Strategy flags -------------------------------------------
        recursionDepth := 1,
        recursionMaxDepth := 5,
        maximalBoundOnTopDegree := 50,
        useMemory := true,                        # whether to use memory for group elements
        # -- Action oracles -------------------------------------------
        act := OnPoints,
        isEqual := \=,                              # isEqual returns true if p and b are equal
                                                    #         returns false if we cannot determine whether p and b are equal
        isDistinguishable := {a, b} -> true,        # if IsEqual returns false (so a = b could not be proven, but might be true),
                                                    # this returns true if we can distinguish between a and b.
                                                    # If points are not distinguishable, we assume they are equal for the moment
        isConflict := {p, indistinguishablePoints} -> false,     # return true if domain might not be invariant under the action
        resolveConflict := {p, indistinguishablePoints} -> fail, # try to find different domain or return fail,
    );

    # ==================================================
    # Input validation
    # ==================================================

    # Check recog node
    if not IsRecogNode(ri) then
        Error("First argument must be a recog node");
    fi;

    # Validate and apply user-provided options
    for name in RecNames(userOptions) do
        if not name in RecNames(options) then
            Error("Invalid option name: ", name);
        fi;
        options.(name) := userOptions.(name);
    od;

    # ==================================================
    # Algorithm
    # ==================================================

    if options.recursionDepth > options.recursionMaxDepth then
        return fail;
    fi;
    Info(WPR_Info, 2, "Recursion Depth: ", options.recursionDepth);
    if options.useMemory then
        G := Group(ri!.gensHmem);
    else
        G := Grp(ri);
    fi;
    t := [One(G)];
    domain := [s];
    i := 1;

    # -- Iteration -------------------------------------
    while i <= Length(domain) do
        a := domain[i];
        for g in GeneratorsOfGroup(G) do
            p := options.act(a, StripMemory(g));
            foundNew := true;
            indistinguishablePoints := [];
            for j in [1 .. Length(domain)] do
                b := domain[j];
                isEq := options.isEqual(p, b);
                if isEq then
                    foundNew := false;
                    break;
                # can we distinguish between p and b via our action?
                elif not options.isDistinguishable(p, b) then
                    Add(indistinguishablePoints, b);
                    # We cannot distinguish between these points,
                    # so we need to assume for now that they are equal.
                    foundNew := false;
                    # If we can prove a conflict, we need to adjust the domain
                    if options.isConflict(p, indistinguishablePoints) then
                        Info(WPR_Info, 2, "Found conflict");
                        # The domain is not invariant under the action,
                        # thus we try compute a new point to resolve the conflict
                        s := options.resolveConflict(p, indistinguishablePoints);
                        if s = fail then
                            Info(WPR_Info, 2, "Could not resolve conflict");
                            return fail;
                        fi;
                        Error("Examine Me!");
                        options.recursionDepth := options.recursionDepth + 1;
                        return BlackBoxOrbit(ri, s, options);
                    fi;
                fi;
            od;
            if foundNew then
                if Length(domain) >= options.maximalBoundOnTopDegree then
                    Info(WPR_Info, 2, "Domain too large, abort now");
                    return fail;
                fi;
                Add(t, t[i] * g);
                Add(domain, p);
            fi;
        od;
        i := i + 1;
    od;
    return rec(domain := domain, transversal := t);
end);