#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  TopGroupDomainAndAction/ViaLinearAction.gi
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
##  This file contains functions to compute a top group domain
##  and the induced action considering subgroups in product action
##  and their orbits under 1.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_TopGroupDomainViaProductAction", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupDomainViaProductAction(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> Set(Orbit(a, 1)) = Set(Orbit(b, 1));
    userOptions.isDistinguishable := {a, b} -> Size(Intersection(Orbit(a, 1), Orbit(b, 1))) = 1;
    userOptions.isConflict := {p, P} -> Length(P) >= 1;
    userOptions.resolveConflict := {p, P} -> Group(Concatenation(List(Concatenation(P, [p]), GeneratorsOfGroup)));

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("WPR_TopGroupActionViaProductAction", function(args...)
    local g, domain, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupActionViaProductAction(g, domain[, options])");
    fi;

    g := args[1];
    domain := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> Set(Orbit(a, 1)) = Set(Orbit(b, 1));

    return BlackBoxAction(g, domain, userOptions);
end);