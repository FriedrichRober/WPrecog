BindGlobal("WPR_InitOptions", function(ri, data, options)
    local timer, name, G, isPrim, N, F, M, d;
    G := Grp(ri);
	Info(WPR_Info, 1, "Step 0: Initialize bounds for wreath product recognition");
	Info(WPR_Info, 2, "--------------------------------------------------------");
	if IsPermGroup(G) then
        isPrim := fail;
        if options.action = fail and options.checkPrimitivity then
            isPrim := IsPrimitive(G);
        elif options.action = "product action" then
            isPrim := true;
        elif options.action = "imprimitive action" then
            isPrim := false;
        fi;
        if isPrim = true or isPrim = fail then
            # assume product action
            N := NrMovedPoints(G);
            F := Collected(PartialFactorization(N, 0));
            M := F[1][2]; # upper bound for top degree
            if ForAll(F, f -> f[2] = M) then
                options.action := "product action";
                isPrim := true;
            fi;
        fi;
        if isPrim = false or isPrim = fail then
            N := NrMovedPoints(G);
            M := N / First(Primes, p -> RemInt(N, p) = 0); # upper bound for top degree
            options.action := "imprimitive action";
            isPrim := false;
        fi;
	elif IsMatrixGroup(G) then
        d := DimensionOfMatrixGroup(G);
        M := d / First(Primes, p -> RemInt(d, p) = 0); # upper bound for top degree
        options.action := "imprimitive action";
    else
        Error("todo");
	fi;

    options.ForSingleComponentGroup.M := Minimum(options.ForSingleComponentGroup.M, M);
    return true;
end);

BindGlobal("WPR_SingleComponentGroup", function(ri, data, options)
    local timer, name, S;
    # ------------------------------------------------------------- #
	# Step 1  : Compute a single-component group S isomorphic to T  #
	# ------------------------------------------------------------- #
	Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step 1 : ", "Compute a single-component group S");
	Info(WPR_Info, 2, "---------------------------------------------------------------");
    timer := Runtime();
	S := SingleComponentGroup(ri, options.ForSingleComponentGroup); # with memory
    timer := Runtime() - timer;
    Info(WPR_Info, 1, "Time : ", timer / 1000.0, " seconds");
    data.S := S;
    return true;
end);

BindGlobal("WPR_TopGroupDomain", function(ri, data, options)
    local timer, name, G, S, S_withoutMem, domainData, W;
    G := Grp(ri);
    S := data.S;
    S_withoutMem := Group(List(GeneratorsOfGroup(S), StripMemory));
    Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step 3 : ", "Compute a faithful H-set");
	Info(WPR_Info, 2, "-------------------------------------------------------------------");
    timer := Runtime();
    # exploit representations
    if IsPermGroup(G) then
        if options.action = "product action" then
	        domainData := TopGroupDomainViaProductAction(ri, S_withoutMem, options.ForTopGroupDomain);
        elif options.action = "imprimitive action" then
            # TODO
        fi;
    elif IsMatrixGroup(G) then
        if options.action = "imprimitive action" then
            W := TopGroupPointViaLinearAction(S);
            data.invariantSubspace := W;
            domainData := TopGroupDomainViaLinearAction(ri, W, options.ForTopGroupDomain);
        fi;
    fi;
    # generic fallback
    if not IsBound(domainData) then
        if options.simpleBaseComp then
            domainData := TopGroupDomainViaSimpleSubgroup(ri, S_withoutMem, options.ForTopGroupDomain);
        else
            domainData := TopGroupDomainViaGenericSubgroup(ri, S_withoutMem, options.ForTopGroupDomain);
        fi;
    fi;
    timer := Runtime() - timer;
    if domainData = fail then
        Info(WPR_Info, 1, "Failure: invalid domain");
        Info(WPR_Info, 1, "Time : ", timer / 1000.0, " seconds");
        return fail;
    fi;
    for name in RecNames(domainData) do
        data.(name) := domainData.(name);
        data.m := Length(domainData.domain);
    od;
	Info(WPR_Info, 2, "m = ", data.m);
    Info(WPR_Info, 1, "Time : ", timer / 1000.0, " seconds");
    return true;
end);