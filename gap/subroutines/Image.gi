#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  Image.gi
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
##  This file contains all high-level functions related to the embedding
##  and image computations. It delegates the actual workload to functions
##  contained in subroutines/Image based on the input group and the options.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_SetupImage", function(ri, data, options)
    local G;
    G := Grp(ri);
    # exploit representations
    if IsPermGroup(G) then
        if options.action = "product action" then
	        return SetupImageViaProductAction(ri, data, options);
        elif options.action = "imprimitive action" then
            return SetupImageViaImprimitiveAction(ri, data, options);
        fi;
    elif IsMatrixGroup(G) then
        if options.action = "imprimitive action" then
            return SetupImageViaLinearAction(ri, data, options);
        fi;
    fi;
    # generic fallback
    return SetupImageViaConjugationAction(ri, data, options);
end);


BindGlobal("WPR_Image", function(g, ri, data, options)
    local G;
    G := Grp(ri);
    # exploit representations
    if IsPermGroup(G) then
        if options.action = "product action" then
	        return ImageViaPermutationAction(g, ri, data, options);
        elif options.action = "imprimitive action" then
            return ImageViaPermutationAction(g, ri, data, options);
        fi;
    elif IsMatrixGroup(G) then
        if options.action = "imprimitive action" then
            return ImageViaLinearAction(g, ri, data, options);
        fi;
    fi;
    # generic fallback
    return ImageViaConjugationAction(g, ri, data, options);
end);

BindGlobal("WPR_Projection", function(g, ri, data, options)
    local G;
    G := Grp(ri);
    # exploit representations
    if IsPermGroup(G) then
        if options.action = "product action" then
	        return WPR_TopGroupActionViaProductAction(g, data.domain);
        elif options.action = "imprimitive action" then
            # TODO
        fi;
    elif IsMatrixGroup(G) then
        if options.action = "imprimitive action" then
            return WPR_TopGroupActionViaLinearAction(g, data.domain);
        fi;
    fi;
    # generic fallback
    return TopGroupActionViaSubgroup(g, data.domain);
end);

BindGlobal("WPR_Embedding", function(ri, data, options)
    local res, W, embedding;
    Info(WPR_Info, 2, "\n");
	Info(WPR_Info, 1, "Step ", data.currentStep, ": Compute embedding function");
	Info(WPR_Info, 2, "--------------------------------------------------");
    data.currentStep := data.currentStep + 1;
	data.projectionFunction  := function(g)
		return WPR_Projection(StripMemory(g), ri, data, options);
	end;

    # Error("Halt!");
    res := WPR_SetupImage(ri, data, options);
    if res = fail then
        return fail;
    fi;
    data.imageFunc := function(g)
		return WPR_Image(StripMemory(g), ri, data, options);
	end;

    W := data.parentWreathProduct;
    embedding := MappingByFunction(Grp(ri), W, data.imageFunc);
    if embedding = fail then
        return fail;
    fi;

    data.embedding := embedding;

    return true;
end);
