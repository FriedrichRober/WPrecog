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

BindGlobal( "RecogniseWreathProduct", function(args...)
    local ri, userOptions, options, name, timer, data, res, output;

    # =======================================================================
    # Input extraction
    # =======================================================================

    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: RecogniseWreathProduct(ri[, options])");
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
        forSingleComponentGroup := rec(),
        forTopGroupDomain := rec(),
        action := fail,
        checkPrimitivity := false,
        assumeSimpleBaseComponent := true,
        recogniseComponents := true,
        recogniseBaseComponentBeforeDomain := true,
        recogniseBaseComponentViaIsomorphism := false,
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
    res := WPR_InitOptions(ri, data, options);
    if res = fail then
        return output;
    fi;

    # Single Component Group
    res := WPR_SingleComponentGroup(ri, data, options);
    if res = fail then
        return output;
    fi;

    # Recognise Base Component?
    if options.recogniseComponents
        and options.recogniseBaseComponentBeforeDomain
    then
        WPR_RecogniseBaseComponent(ri, data, options);
    fi;

    # Top Group Domain
	res := WPR_TopGroupDomain(ri, data, options);
    if res = fail then
        return output;
    fi;

    # TODO: handle this cleaner
    if IsMatrixGroup(Grp(ri)) then
        if options.action = "imprimitive action" then
            if Dimension(data.domain[1]) * data.m = DimensionOfMatrixGroup(Grp(ri)) then
                data.B := Concatenation(List(data.domain, Basis))^(-1);
                output.res := true;
                return output;
            fi;
        fi;
        Error("TODO");
    fi;

    # Recognise Base Component?
    if options.recogniseComponents
        and not options.recogniseBaseComponentBeforeDomain
    then
        WPR_RecogniseBaseComponent(ri, data, options);
    fi;

    # Embedding
    res := WPR_Embedding(ri, data, options);
    if res = fail then
        return output;
    fi;

    timer := Runtime() - timer;
    Info(WPR_Info, 1, "===============================================================");
    Info(WPR_Info, 1, "Total Time: ", FormatFloat(timer / 1000.0), " seconds");

    output.res := true;
    return output;
end);
