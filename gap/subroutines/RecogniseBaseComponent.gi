BindGlobal("WPR_RecogniseBaseComponent", function(ri, data, options)
    local timer, name, G, gensS, S_withoutMem, isoReductionForS, S_reduced, riS_reduced, orb;
    G := Grp(ri);
    Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step ", data.currentStep, ": Recognise base component");
	Info(WPR_Info, 2, "-------------------------------------------------------------------");
    data.currentStep := data.currentStep + 1;
    timer := Runtime();
    gensS := List(GeneratorsOfGroup(data.S), StripMemory);
    S_withoutMem := Group(gensS);
    isoReductionForS := GroupHomomorphismByFunction(S_withoutMem, S_withoutMem, x -> x);
    if IsPermGroup(G) then
        if options.action = "product action" then
            orb := Orbit(S_withoutMem, 1);
            isoReductionForS := ActionHomomorphism(S_withoutMem, orb, OnPoints);
        fi;
    fi;
    data.isoReductionForS := isoReductionForS;
    S_reduced := Image(data.isoReductionForS);
    if options.recogniseBaseComponentViaIsomorphism then
        Error("TODO");
    else
        riS_reduced := RecogniseGroup(S_reduced);
    fi;
    timer := Runtime() - timer;
    if riS_reduced = fail then
        Info(WPR_Info, 1, "Failure: could not recognise base component");
        Info(WPR_Info, 1, "Time : ", timer / 1000.0, " seconds");
        return fail;
    fi;
    data.riS_reduced := riS_reduced;
    data.niceGensForS := ResultOfStraightLineProgram(slptonice(data.riS_reduced), gensS);
    Info(WPR_Info, 1, "Time : ", timer / 1000.0, " seconds");
    return true;
end);