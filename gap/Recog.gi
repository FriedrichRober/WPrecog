
BindGlobal( "RecogniseWreathProduct", function(args...)
    local ri, userOptions, options, name, data, res, output;

    # ==================================================
    # Input extraction
    # ==================================================
    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: RecogniseWreathProduct(ri[, options])");
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
        ForSingleComponentGroup := rec(),
        ForTopGroupDomain := rec(),
        action := fail,
        checkPrimitivity := false,
        simpleBaseComp := true,
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
    data := rec();
    output := rec(res := fail, data := data, options := options);

    # Init options
    res := WPR_InitOptions(ri, data, options);
    if res = fail then
        return output;
    fi;

    # Single Component Group
    res := WPR_SingleComponentGroup(ri, data, options);
    if res = fail then
        return output;
    fi;

    # Top Group Domain
	res := WPR_TopGroupDomain(ri, data, options);
    if res = fail then
        return output;
    fi;

    output.res := true;
    return output;
end);
