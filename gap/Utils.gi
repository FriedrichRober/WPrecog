#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  Utils.gi
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
##  This file contains utility functions.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################


FormatFloat := function(args...)
    local x, nrDigits, s, parts, decimals;
    if Length(args) = 0 or Length(args) > 2 then
        Error("Usage: SingleComponentGroup(x[, nrDigits])");
    fi;

    x := args[1];
    nrDigits := 3;
    if Length(args) = 2 then
        nrDigits := args[2];
    fi;

    s := String(x);
    parts := SplitString(s, ".");
    if Length(parts) = 1 then
        return Concatenation(s, ".000");
    fi;
    decimals := parts[2];
    while Length(decimals) < nrDigits do
        decimals := Concatenation(decimals, "0");
    od;
    if Length(decimals) > nrDigits then
        decimals := decimals{[1..nrDigits]};
    fi;
    return Concatenation(parts[1], ".", decimals);
end;
