#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  SingleComponentGroup.gi
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
##  This file contains the functions to blindly descend into
##  a single-component group.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################


#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Called by WreathProductDecomposition.
## Wrapper for SingleComponentGroup.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_SingleComponentGroup", function(ri, data, options)
    local timer, name, singleCompData;
	Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step ", data.currentStep, ": Compute a single-component group S");
	Info(WPR_Info, 2, "---------------------------------------------------------------");
    data.currentStep := data.currentStep + 1;
    timer := Runtime();
	singleCompData := SingleComponentGroup(ri, options.forSingleComponentGroup); # with memory
    timer := Runtime() - timer;
    Info(WPR_Info, 1, "Time: ", FormatFloat(timer / 1000.0), " seconds");
    for name in RecNames(singleCompData) do
        data.(name) := singleCompData.(name);
    od;
    return true;
end);

#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Args: ri[, options]
## - ri is recog node
## - options are user options, see in code
##
## Returns a supposedly single-component group.
## The backbone of WPR_SingleComponentGroup.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

InstallGlobalFunction("SingleComponentGroup", function(args...)
    local ri, userOptions, options, name, S, l1, l2, primes, R, B, i, descendOptions, logEps;

    # =======================================================================
    # Input extraction
    # =======================================================================

    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: SingleComponentGroup(ri[, options])");
    fi;

    ri := args[1];
    userOptions := rec();
    if Length(args) = 2 then
        userOptions := args[2];
    fi;

    # =======================================================================
    # Default options
    # =======================================================================
    options := rec(
        # -- Error control -------------------------------------------------------------------------
        heuristic := true,                      # whether to use heuristic values instead of explicit ones
        lazy := true,                           # whether to try at most 2 log(m) tries
        eps := 1/100,                           # error probability
        M := 20,                                # upper bound for degree of top group

        # -- Strategy flags ------------------------------------------------------------------------
        primes := [2, 3, 5, 7],                 # prime numbers to use in a step
        randomElementGenerator                  # function to generate random elements
          := QuickRandomSubproduct,
        useMemory := true,                      # whether to use memory for group elements
        useBaseGroupInPhase2 := false,          # use base group instead of previous group in phase 2
        useBaseGroupAtEnd := true,              # use base group at end (ignored if useBaseGroupInPhase2)
        socleType := fail,                      # socleType
        updateParameters := false,              # whether to update parameters after phase 1

        # -- Phase iteration counts ----------------------------------------------------------------
        l1 := 5,                                # number of iterations in phase 1
        l2 := 35,                               # number of iterations in phase 2
        l1Min := 5,                             # minimal number of iterations in phase 1
        l2Min := 5,                             # minimal number of iterations in phase 2

        # -- Normal closure effort parameters ------------------------------------------------------
        n1 := 2,                                # effort per step in phase 1
        n2 := 2,                                # effort per step in phase 2
        n0 := 2,                                # effort at final step (used only if useBaseGroupAtEnd)
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

    if options.useMemory then
        S := Group(ri!.gensHmem);
    else
        S := Grp(ri);
    fi;
    R := WPR_SCG_ReductionRate(options.socleType);
    if not options.heuristic then
        primes := [2];
        Info(WPR_Info, 2, "Reduction rate probabilities for a step in each phase:");
        Info(WPR_Info, 2, "\tP_1 = ", R[1]);
        Info(WPR_Info, 2, "\tP_2 = ", R[2]);
    else
        primes := options.primes;
    fi;

    # -- Phase 1 ------------------------------------------------------------
    l1 := WPR_SCG_NrIterations1(R[1], options.eps/2);
    if options.heuristic then
        l1 := Minimum(l1, options.l1);
        l1 := Maximum(l1, options.l1Min);
    fi;
    descendOptions := rec(
        randomElementGenerator := options.randomElementGenerator,
        n := options.n1,
        primes := primes,
    );
    Info(WPR_Info, 2, "Number of iterations for Phase 1: ", l1);
	Info(WPR_Info, 3, "Start iteration for Phase 1...");
    for i in [1 .. l1] do
        S := WPR_SCG_Descend(S, S, descendOptions);
    od;
    Info(WPR_Info, 3, "...Finished iteration for Phase 1");

    # -- Update Parameters --------------------------------------------------
    B := S; # supposedly inside base group
    if options.updateParameters then
        WPR_SCG_UpdateParameters(options, B);
    fi;

    # -- Phase 2 ------------------------------------------------------------
    l2 := WPR_SCG_NrIterations2(R[2], options.M, options.eps/2);
    if options.lazy then
        l2 := Minimum(l2, Int(Ceil(2*Log(Float(options.M)))));
    fi;
    if options.heuristic then
        l2 := Minimum(l2, options.l2);
        l2 := Maximum(l2, options.l2Min);
    fi;
    descendOptions := rec(
        randomElementGenerator := options.randomElementGenerator,
        n := options.n2,
        primes := primes,
    );
    Info(WPR_Info, 2, "Number of iterations for Phase 2: ", l2);
    Info(WPR_Info, 3, "Start iteration for Phase 2...");
    for i in [1 .. l2] do
        if options.useBaseGroupInPhase2 then
            S := WPR_SCG_Descend(S, B, descendOptions);
        else
            S := WPR_SCG_Descend(S, S, descendOptions);
        fi;
    od;
    Info(WPR_Info, 3, "...Finished iteration for Phase 2");

    # -- Post ---------------------------------------------------------------
    # As a pre-caution, take normal closure in base group, if we descended too far
    if not options.useBaseGroupInPhase2 and options.useBaseGroupAtEnd then
        Info(WPR_Info, 3, "Start computing normal closure in base group...");
        S := QuickNormalClosure(S, B, options.n0);
        Info(WPR_Info, 3, "... Finished computing normal closure in base group");
    fi;
    return rec(B := B, S := S);
end);
