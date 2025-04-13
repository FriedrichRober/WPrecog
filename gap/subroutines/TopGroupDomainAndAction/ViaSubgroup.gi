#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  TopGroupDomainAndAction/ViaSubgroup.gi
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
##  and the induced action considering subgroups.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("DoCommuteElms", function(x, y)
    return x * y = y * x;
end);

BindGlobal("DoCommuteGroups", function(A, B)
    local a, b;
    for a in GeneratorsOfGroup(A) do
        for b in GeneratorsOfGroup(B) do
            if not DoCommuteElms(a, b) then
                return false;
            fi;
        od;
    od;
    return true;
end);

BindGlobal("WPR_TopGroupDomainViaSimpleSubgroup", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupDomainViaSimpleSubgroup(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> not DoCommuteGroups(a,b); # if groups do not commute, groups must be equal
    userOptions.isDistinguishable := {a, b} -> true; # if groups commute, groups must be disjoint

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("WPR_TopGroupDomainViaGenericSubgroup", function(args...)
    local ri, s, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupDomainViaGenericSubgroup(ri, s[, options])");
    fi;

    ri := args[1];
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> false; # we cannot check for equality without further assumptions
    userOptions.isDistinguishable := {a, b} -> DoCommuteGroups(a,b); # if groups commute, then they might lie in different components
    userOptions.isConflict := {p, P} -> Length(P) >= 2; # we found two groups with which p does not commute
    userOptions.resolveConflict := {p, P} -> Group(Concatenation(List(Concatenation(P, [p]), GeneratorsOfGroup))); # we merge these groups and hope that we now cover the socle in the component

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("TopGroupActionViaSubgroup", function(args...)
    local g, domain, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupActionViaLinearAction(g, domain[, options])");
    fi;

    g := args[1];
    domain := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> not DoCommuteGroups(a,b);

    return BlackBoxAction(g, domain, userOptions);
end);
