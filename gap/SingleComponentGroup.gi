# args: ri[, options]
# - ri is recog node
# - options are user options, see in code
# Returns single-component group
InstallGlobalFunction("SingleComponentGroup", function(args...)
    local ri, userOptions, options, name, S, l1, l2, P, R, B, i, descendOptions, logEps;

    # ==================================================
    # Input extraction
    # ==================================================
    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: SingleComponentGroup(ri[, options])");
    fi;

    ri := args[1];
    userOptions := rec();
    if Length(args) = 2 then
        userOptions := args[2];
    fi;

    # ==================================================
    # Default options
    # ==================================================
    options := rec(
        # -- Error control -------------------------------------------
        heuristic := true,                      # whether to use heuristic values instead of explicit ones
        eps := 1/100,                           # error probability
        M := 20,                                # upper bound for degree of top group

        # -- Strategy flags -------------------------------------------
        P := [2, 3, 5, 7],                      # prime numbers to use in a step
        randomElementGenerator                  # function to generate random elements
          := QuickRandomSubproduct,
        useMemory := true,                     # whether to use memory for group elements
        useBaseGroupInPhase2 := false,          # use base group instead of previous group in phase 2
        useBaseGroupAtEnd := true,              # use base group at end (ignored if useBaseGroupInPhase2)
        socleType := fail,                      # socleType
        updateParameters := false,              # whether to update parameters after phase 1

        # -- Phase iteration counts -------------------------------------------
        l1 := 5,                                # number of iterations in phase 1
        l2 := 35,                               # number of iterations in phase 2

        # -- Normal closure effort parameters ---------------------------------
        n1 := 5,                                # effort per step in phase 1
        n2 := 5,                                # effort per step in phase 2
        n0 := 2,                                # effort at final step (used only if useBaseGroupAtEnd)
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

    if options.useMemory then
        S := Group(ri!.gensHmem);
    else
        S := Grp(ri);
    fi;
    R := SingleComponentGroup_ReductionRate(options.socleType);
    if not options.heuristic then
        P := [2];
        Info(WPR_Info, 2, "Reduction rate probabilities for a step in each phase:");
        Info(WPR_Info, 2, "\tP_1 = ", R[1]);
        Info(WPR_Info, 2, "\tP_2 = ", R[2]);
    else
        P := options.P;
    fi;

    # -- Phase 1 -------------------------------------------
    l1 := SingleComponentGroup_NrIterations1(R[1], options.eps/2);
    if options.heuristic then
        l1 := Minimum(l1, options.l1);
    fi;
    descendOptions := rec(
        randomElementGenerator := options.randomElementGenerator,
        n := options.n1,
        P := P,
    );
    Info(WPR_Info, 2, "Bounds for Phase 1:");
	Info(WPR_Info, 2, "\tl_1 = ", l1);
	Info(WPR_Info, 3, "Start Iteration for Phase 1...");
    for i in [1 .. l1] do
        S := SingleComponentGroup_Descend(S, S, descendOptions);
    od;
    Info(WPR_Info, 3, "...Finished Iteration for Phase 1");

    # -- Update Parameters ---------------------------------
    B := S; # supposedly inside base group
    if options.updateParameters then
        SingleComponentGroup_UpdateParameters(options, B);
    fi;

    # -- Phase 2 -------------------------------------------
    l2 := SingleComponentGroup_NrIterations2(R[2], options.M, options.eps/2);
    if options.heuristic then
        l2 := Minimum(l2, Int(Ceil(2*Log(Float(options.M)))), options.l2);
    fi;
    descendOptions := rec(
        randomElementGenerator := options.randomElementGenerator,
        n := options.n2,
        P := P,
    );
    Info(WPR_Info, 2, "Bounds for Phase 2:");
	Info(WPR_Info, 2, "\tl_2 = ", l2);
    Info(WPR_Info, 3, "Start Iteration for Phase 2...");
    for i in [1 .. l2] do
        if options.useBaseGroupInPhase2 then
            S := SingleComponentGroup_Descend(S, B, descendOptions);
        else
            S := SingleComponentGroup_Descend(S, S, descendOptions);
        fi;
    od;
    Info(WPR_Info, 3, "...Finished Iteration for Phase 2");

    # -- Post -------------------------------------------
    # As a pre-caution, take normal closure in base group, if we descended too far
    if not options.useBaseGroupInPhase2 and options.useBaseGroupAtEnd then
        Info(WPR_Info, 3, "Start computing normal closure in base group...");
        S := QuickNormalClosure(S, B, options.n0);
        Info(WPR_Info, 3, "... Finished computing normal closure in base group");
    fi;
    return S;

end);

# args: H, G, options
# - for options, see SingleComponentGroup.
# - H is current group in iteration step
# - computes normal closure in G
# Returns next group in iteration step
InstallGlobalFunction("SingleComponentGroup_Descend", function(H, G, options)
    local y, ord, P, p, z;
    y := options.randomElementGenerator(H);
    ord := Order(y);
    P := options.P{QuickRandomPermList(Length(options.P))};
    for p in P do
        if (p = 2 and IsEvenInt(ord))
          or RemInt(ord, p) = 0 then
            z := y ^ (ord / p);
            return QuickNormalClosure([z], G, options.n);
        fi;
    od;
    return H;
end);

InstallGlobalFunction( "SingleComponentGroup_ReductionRate", function(socleType)
    if socleType = fail then
        socleType := "Alt";
    fi;
	if socleType = "Alt" then
		return [1/2, 1/3];
	fi;
	Error("Invalid socle type: ", socleType);
end);

InstallGlobalFunction( "SingleComponentGroup_NrIterations1", function(p1, eps)
    local logEps;
    logEps := Log(Float(1 / eps));
    return Int(Ceil(logEps/Log(Float(1 / (1 - p1)))));
end);

InstallGlobalFunction( "SingleComponentGroup_NrIterations2", function(p2, M, eps)
    local delta, c, r, k, m, logDelta, l_lin, l_log, l2;
    delta := eps/2;
    logDelta := Log(Float(1 / delta));
    c := 2 * (1 + 1/p2);
    r := 1 - p2/2;
    m := 2*c;
    if m < M then
        k := Int(Ceil(Log(Float(M/c))/Log(Float(1/r))));
        l_log := Int(Ceil(Maximum(Float(4 * k), 16 * logDelta)));
    else
        l_log := 0;
    fi;
    l_lin := Int(Ceil(Maximum(Float(2 / p2 * m), 8 / p2 * logDelta)));
    l2 := Int(Ceil(Maximum(Float(2 / p2 * M), 8 / p2 * logDelta)));
    return Minimum(l2, l_log + l_lin);
end);

InstallGlobalFunction("SingleComponentGroup_UpdateParameters", function(options, B)
    return; # TODO
end);

