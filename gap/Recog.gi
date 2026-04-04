#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  Recog.gi
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
##  This file contains the main function for recognising wreath products.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal( "WPR_ConvertOutput", function(output, debug)
    if debug then
        return output;
    elif output.res = true then
        return output.data.embedding;
    else
        return fail;
    fi;
end);

BindGlobal( "WreathProductDecomposition", function(args...)
    local ri, userOptions, options, name, timer, data, res, output;

    # =======================================================================
    # Input extraction
    # =======================================================================

    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: WreathProductDecomposition(ri[, options])");
    fi;

    if IsGroup(args[1]) then
        ri := RecogNode(args[1]);
    elif IsRecogNode(args[1]) then
        ri := args[1];
    else
        Error("Usage: first argument must be a group or recog node");
    fi;
    userOptions := rec();
    if Length(args) = 2 then
        userOptions := args[2];
    fi;

    # =======================================================================
    # Default options
    # =======================================================================

    options := rec(
        forSingleComponentGroup := rec(),
        forTopGroupDomain := rec(),
        action := fail,
        checkPrimitivity := false,
        assumeSimpleBaseComponent := true,
        recogniseBaseComponent := true,
        recogniseBaseComponentBeforeDomain := true,
        recogniseBaseComponentViaIsomorphism := false,
        M := fail,
        debug := false,
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

    data := rec(currentStep := 0);
    output := rec(res := fail, data := data, options := options);
    timer := Runtime();

    # Init options
    output.res := WPR_InitOptions(ri, data, options);
    if output.res = fail then
        return WPR_ConvertOutput(output, options.debug);
    fi;

    # Single Component Group
    output.res := WPR_SingleComponentGroup(ri, data, options);
    if output.res = fail then
        return WPR_ConvertOutput(output, options.debug);
    fi;

    # Recognise Base Component?
    if options.recogniseBaseComponent
        and options.recogniseBaseComponentBeforeDomain
    then
        WPR_RecogniseBaseComponent(ri, data, options);
    fi;

    # Top Group Domain
	output.res := WPR_TopGroupDomain(ri, data, options);
    if output.res = fail then
        return WPR_ConvertOutput(output, options.debug);
    fi;

    # TODO: handle this cleaner
    if IsMatrixGroup(Grp(ri)) then
        if options.action = "imprimitive action" then
            if Dimension(data.domain[1]) * data.m = DimensionOfMatrixGroup(Grp(ri)) then
                data.B := Concatenation(List(data.domain, Basis))^(-1);
                output.res := true;
                return WPR_ConvertOutput(output, options.debug);
            fi;
        fi;
        Error("TODO");
    fi;

    # Recognise Base Component?
    if options.recogniseBaseComponent
        and not options.recogniseBaseComponentBeforeDomain
    then
        WPR_RecogniseBaseComponent(ri, data, options);
    fi;

    # Embedding
    output.res := WPR_Embedding(ri, data, options);
    if output.res = fail then
        return WPR_ConvertOutput(output, options.debug);
    fi;

    timer := Runtime() - timer;
    Info(WPR_Info, 1, "===============================================================");
    Info(WPR_Info, 1, "Total Time: ", FormatFloat(timer / 1000.0), " seconds");

    output.res := true;
    return WPR_ConvertOutput(output, options.debug);
end);
