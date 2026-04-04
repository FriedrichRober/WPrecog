#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  TopGroupDomainAndAction.gi
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##
##  This file's authors include Friedrich Rober.
##
##  Please refer to the COPYRIGHT file for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##
##  This file contains the black box functions to compute
##  a top group domain and the induced action.
##  Specialized functions are in the respective directory.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################


#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Called by WreathProductDecomposition.
## Wrapper for WPR_TopGroupDomainVia___.
## These use BlackBoxOrbit under the hood.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_TopGroupDomain", function(ri, data, options)
    local timer, name, G, gensS, S_withoutMem, domainData, W;
    G := Grp(ri);
    gensS := List(GeneratorsOfGroup(data.S), StripMemory);
    if IsBound(data.niceGensForS) and Length(data.niceGensForS) < Length(gensS) then
        gensS := List(data.niceGensForS, StripMemory);
    fi;
    S_withoutMem := Group(gensS);
    Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step ", data.currentStep, ": Compute a faithful H-set");
	Info(WPR_Info, 2, "-------------------------------------------------------------------");
    data.currentStep := data.currentStep + 1;
    timer := Runtime();
    # exploit representations
    if IsPermGroup(G) then
        if options.action = "product action" then
	        domainData := WPR_TopGroupDomainViaProductAction(ri, S_withoutMem, options.forTopGroupDomain);
        elif options.action = "imprimitive action" then
            # TODO
        fi;
    elif IsMatrixGroup(G) then
        if options.action = "imprimitive action" then
            W := WPR_TopGroupPointViaLinearAction(S_withoutMem);
            data.invariantSubspace := W;
            domainData := WPR_TopGroupDomainViaLinearAction(ri, W, options.forTopGroupDomain);
        fi;
    fi;
    # generic fallback
    if not IsBound(domainData) then
        if options.assumeSimpleBaseComponent then
            domainData := WPR_TopGroupDomainViaSimpleSubgroup(ri, S_withoutMem, options.forTopGroupDomain);
        else
            domainData := WPR_TopGroupDomainViaGenericSubgroup(ri, S_withoutMem, options.forTopGroupDomain);
        fi;
    fi;
    timer := Runtime() - timer;
    if domainData = fail then
        Info(WPR_Info, 1, "Failure: invalid domain");
        Info(WPR_Info, 1, "Time: ", FormatFloat(timer / 1000.0), " seconds");
        return fail;
    fi;
    for name in RecNames(domainData) do
        data.(name) := domainData.(name);
        data.m := Length(domainData.domain);
    od;
	Info(WPR_Info, 2, "m = ", data.m);
    Info(WPR_Info, 1, "Time: ", FormatFloat(timer / 1000.0), " seconds");
    return true;
end);

#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Args: ri, s[, options]
## - ri is recog node with G = Grp(ri)
## - s is a single point on which G acts on
## - options are user options, see in code
##
## Returns a supposedly orbit of s under G.
## The backbone of TopGroupDomainVia___.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

InstallGlobalFunction("BlackBoxOrbit", function(args...)
    local ri, s, userOptions, options, name, G, t, domain, i, a, b, g, p, isEq, foundNew, indistinguishablePoints, j;

    # =======================================================================
    # Input extraction
    # =======================================================================
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: BlackBoxOrbit(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;

    # =======================================================================
    # Default options
    # =======================================================================
    options := rec(
        # -- Strategy flags -------------------------------------------
        recursionDepth := 1,                        # current recursion depth
        recursionMaxDepth := 5,                     # maximal recursion depth
        maximalBoundOnTopDegree := 50,              # if the computed domain exceeds this bound, then return fail
        useMemory := true,                          # whether to use memory for group elements
        # -- Action oracles -------------------------------------------
        act := OnPoints,
        isEqual := \=,                              # isEqual returns true if p and b are equal
                                                    #         returns false if we cannot determine whether p and b are equal
        isDistinguishable := {a, b} -> true,        # Called after IsEqual returns false.
                                                    # We use a heuristic and our action to distinguish between a and b.
                                                    # This returns true if we can distinguish between a and b.
                                                    # If points are not distinguishable, we assume they are equal for the moment
        isConflict := {p, indistinguishablePoints} -> false,     # return true if domain might not be invariant under the action
        resolveConflict := {p, indistinguishablePoints} -> fail, # try to find different domain or return fail,
    );

    # =======================================================================
    # Input validation
    # =======================================================================

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

    # =======================================================================
    # Algorithm
    # =======================================================================

    if options.recursionDepth > options.recursionMaxDepth then
        return fail;
    fi;
    Info(WPR_Info, 3, "Recursion Depth: ", options.recursionDepth);
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
                        Info(WPR_Info, 3, "Found conflict");
                        # The domain is not invariant under the action,
                        # thus we try compute a new point to resolve the conflict
                        s := options.resolveConflict(p, indistinguishablePoints);
                        if s = fail then
                            Info(WPR_Info, 3, "Could not resolve conflict");
                            return fail;
                        fi;
                        options.recursionDepth := options.recursionDepth + 1;
                        return BlackBoxOrbit(ri, s, options);
                    fi;
                fi;
            od;
            if foundNew then
                if Length(domain) >= options.maximalBoundOnTopDegree then
                    Info(WPR_Info, 3, "Domain too large, abort now");
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

#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Args: g, domain[, options]
## - g is a group element
## - domain is a supposedly action domain of g
## - options are user options, see in code
##
## Returns a supposedly action of g on domain.
## The backbone of TopGroupActionVia___.
## Similar to BlackBoxOrbit.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

InstallGlobalFunction("BlackBoxAction", function(args...)
    local g, domain, s, userOptions, options, name, m, images, i, j, a, b, p;

    # =======================================================================
    # Input extraction
    # =======================================================================

    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: BlackBoxAction(g, domain[, options])");
    fi;

    g := args[1];
    domain := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;

    # =======================================================================
    # Default options
    # =======================================================================

    options := rec(
        # -- Action oracles -------------------------------------------
        act := OnPoints,
        isEqual := \=,                              # isEqual returns true if p and b are assumed to be equal
    );

    # =======================================================================
    # Input validation
    # =======================================================================

    # Validate and apply user-provided options
    for name in RecNames(userOptions) do
        if not name in RecNames(options) then
            Error("Invalid option name: ", name);
        fi;
        options.(name) := userOptions.(name);
    od;

    # =======================================================================
    # Algorithm
    # =======================================================================

    m := Length(domain);
	images := EmptyPlist(m);
    for i in [1 .. m] do
        a := domain[i];
        p := options.act(a, g);
        for j in [1 .. m] do
            if j in images then
                continue;
            fi;
            b := domain[j];
            if options.isEqual(p, b) then
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
