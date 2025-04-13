#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  SingleComponentGroup/Helpers.gi
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
##  This file contains helper functions for SingleComponentGroup.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################


#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Args: H, G, options
## - H is current group in iteration step
## - computes normal closure in G
## - for options, see SingleComponentGroup.
##
## Returns next group in iteration step.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_SCG_Descend", function(H, G, options)
    local y, ord, primes, p, z;
    y := options.randomElementGenerator(H);
    ord := Order(y);
    primes := options.primes{QuickRandomPermutedList(Length(options.primes))};
    for p in primes do
        if (p = 2 and IsEvenInt(ord))
          or RemInt(ord, p) = 0 then
            z := y ^ (ord / p);
            return QuickNormalClosure([z], G, options.n);
        fi;
    od;
    return H;
end);

#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Explicit probability estimations.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal( "WPR_SCG_ReductionRate", function(socleType)
    if socleType = fail then
        socleType := "Alt";
    fi;
	if socleType = "Alt" then
		return [1/2, 1/3];
	fi;
	Error("Invalid socle type: ", socleType);
end);

BindGlobal( "WPR_SCG_NrIterations1", function(p1, eps)
    local logEps;
    logEps := Log(Float(1 / eps));
    return Int(Ceil(logEps/Log(Float(1 / (1 - p1)))));
end);

BindGlobal( "WPR_SCG_NrIterations2", function(p2, M, eps)
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

#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##
## Returns updated parameters after Phase 1,
## assuming that B is a subgroup of the base group.
##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_SCG_UpdateParameters", function(options, B)
    return; # TODO
end);
