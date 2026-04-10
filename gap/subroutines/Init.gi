#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  Init.gi
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
##  This file contains a function to initialize our assumptions
##  and parameters about the wreath product in which our input group
##  naturally fits. It might use hints provided by the User.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_InitOptions", function(ri, data, options)
    local timer, name, G, isPrim, N, F, M, d, q;
    G := Grp(ri);
	Info(WPR_Info, 1, "Step ", data.currentStep, ": Initialize bounds for wreath product recognition");
	Info(WPR_Info, 2, "--------------------------------------------------------");
    data.currentStep := data.currentStep + 1;
	if IsPermGroup(G) then
        Info(WPR_Info, 2, "G is permutation group on ", NrMovedPoints(G), " points");
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
            M := Gcd(List(F, f -> f[2])); # upper bound for top degree
            if options.action = fail and M > 1 and M < 10 then
                options.action := "product action";
                isPrim := true;
            fi;
            if options.action = "product action" then
                Info(WPR_Info, 2, "Assume product action on n^m points with m <= ", M);
            fi;
        fi;
        if isPrim = false or isPrim = fail then
            N := NrMovedPoints(G);
            M := N / First(Primes, p -> RemInt(N, p) = 0); # upper bound for top degree
            options.action := "imprimitive action";
            isPrim := false;
            Info(WPR_Info, 2, "Assume imprimitive action on n*m points with m <= ", M);
        fi;
	elif IsMatrixGroup(G) then
        d := DimensionOfMatrixGroup(G);
        q := Size(FieldOfMatrixGroup(G));
        Info(WPR_Info, 2, "G is matrix group of dimension ", d, " over GF(", q, ")");
        M := d / First(Primes, p -> RemInt(d, p) = 0); # upper bound for top degree
        options.action := "imprimitive action";
        Info(WPR_Info, 2, "Assume imprimitve action on n*m points with m <= ", M);
    else
        Error("todo");
	fi;

    if options.M <> fail then
        Info(WPR_Info, 2, "User specified an upper bound for m");
        options.M := Minimum(M, options.M);
        Info(WPR_Info, 2, "Assume m <= ", M);
    else
        options.M := M;
    fi;
    options.forSingleComponentGroup.M := options.M;
    options.forTopGroupDomain.maximalBoundOnTopDegree := options.M;
    return true;
end);
