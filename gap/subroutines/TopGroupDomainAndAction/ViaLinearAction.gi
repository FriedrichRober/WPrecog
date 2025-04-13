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
##  and the induced action considering subgroups and their linear actions.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPR_TopGroupPointViaLinearAction", function(S)
    local d, q, M, v1, v2, L1, L2, W;
    d := DimensionOfMatrixGroup(S);
    q := Size(DefaultFieldOfMatrixGroup(S));
    # S acts on W + U
    # W invariant space, dim(W) = d
    # U fixed space, dim(U) = (m-1)*d
    M := GModuleByMats(GeneratorsOfGroup(S), GF(q));
    v1 := QuickRandomVector(d,q);
    v2 := QuickRandomVector(d,q);
    # basis of W + < [v]\rho_U >
    L1 := MTX.SubmoduleGModule(M, v1);;
    L2 := MTX.SubmoduleGModule(M, v2);;
    W := Intersection(VectorSpace(GF(q), L1), VectorSpace(GF(q), L2));
    return W;
end);

BindGlobal("WPR_TopGroupDomainViaLinearAction", function(args...)
    local ri, s, q, userOptions;
    if Length(args) < 2 or Length(args) > 3 then
        Error("Usage: WPR_TopGroupDomainViaLinearAction(ri, s[, options])");
    fi;

    ri := args[1];
    q := Size(DefaultFieldOfMatrixGroup(Grp(ri)));
    s := args[2];
    userOptions := rec();
    if Length(args) = 3 then
        userOptions := args[3];
    fi;
    userOptions.isEqual := {a, b} -> a = b;
    userOptions.isDistinguishable := {a, b} -> IsTrivial(Intersection(a, b));
    userOptions.isConflict := {p, P} -> Length(P) >= 1;
    userOptions.resolveConflict := function(p1, P)
        local p2, pNew;
        p2 := P[1];
        pNew := VectorSpace(GF(q), Concatenation(Basis(p1),Basis(p2)));
        Basis(pNew);
        Dimension(pNew);
        return pNew;
    end;

    return BlackBoxOrbit(ri, s, userOptions);
end);

BindGlobal("WPR_TopGroupActionViaLinearAction", function(args...)
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
    userOptions.isEqual := {a, b} -> a = b;

    return BlackBoxAction(g, domain, userOptions);
end);