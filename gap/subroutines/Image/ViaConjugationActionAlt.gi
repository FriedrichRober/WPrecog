#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  ViaConjugationAction.gi
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
##  This file contains functions related to the embedding
##  and image computations. It assumes that the found single-component group
##  is an alternating group which was successfully recognised with
##  standard generators given by (1,2,3) and [ (1,2) ] (3,4,...,n).
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("ImageViaConjugationActionAlt",
function(ri, data, options)
	local t, lambda, n, m, top, base, i, j, k, a, b, aCycles, bCycles, aPoints, bPoints, abPoints, 3Point, 3PosInA, 3PosInB, images;
    t := data.transversal;
    lambda := data.isoS;
	n := NrMovedPoints(Image(lambda));
	m := Length(t);
	top := projFunc(g);
	if top = fail then
		return TemporaryFailure;
	fi;
	base := EmptyPlist(m);
	for i in [1 .. m] do
		j := i ^ top;
		images := EmptyPlist(n);
		a := (stdGensS[1] ^ (t[i] * g * t[j] ^ -1)) ^ lambda;
		if a = fail then
			return TemporaryFailure;
		fi;
		b := (stdGensS[2] ^ (t[i] * g * t[j] ^ -1)) ^ lambda;
		if b = fail then
			return TemporaryFailure;
		fi;
		# aPoints contains images of [1 .. 3] up to cyclic shifting
		aCycles := Cycles(a, [1 .. n]);
		if Set(aCycles, Length) <> [1, 3] then
			return TemporaryFailure;
		fi;
		aPoints := First(aCycles, c -> Length(c) = 3);
		# bPoints contains images of [3 .. n] up to cyclic shifting
		bCycles := Cycles(b, [1 .. n]);
		if IsEvenInt(n) and Set(bCycles, Length) <> [2, n - 2] then
			return TemporaryFailure;
		elif not IsEvenInt(n) and Set(bCycles, Length) <> [1, n -2] then
			return TemporaryFailure;
		fi;
		bPoints := First(bCycles, c -> Length(c) = n - 2);
		# the image of 3 is the unique intersection point of aPoints and bPoints
		abPoints := Intersection(aPoints, bPoints);
		if Length(abPoints) <> 1 then
			return TemporaryFailure;
		fi;
		3Point := abPoints[1];
		images[3] := 3Point;
		3PosInA := Position(aPoints, 3Point);
		images[1] := aPoints[(3PosInA + 1 - 1) mod 3 + 1];
		images[2] := aPoints[(3PosInA + 2 - 1) mod 3 + 1];
		3PosInB := Position(bPoints, 3Point);
		for k in [1 .. n - 3] do
			images[3 + k] := bPoints[(3PosInB + k - 1) mod (n - 2) + 1];
		od;
		base[i] := PermList(images);
	od;
	return Concatenation(base, [top]);
end);