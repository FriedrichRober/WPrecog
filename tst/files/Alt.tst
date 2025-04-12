gap> LoadPackage("WPrecog");
gap> n := 11;;
gap> G := WreathProduct(SymmetricGroup(n), DihedralGroup(IsPermGroup, 8));;
gap> ri := RecogNode(G);;
gap> simpleGroupFamily := "Alt";;
gap> SetInfoLevel(WPR_Info, 3);
gap> WreathProductRecognition(ri, G, simpleGroupFamily);


# Example for slow GAP Socle but fast WPrecog Socle (Factor)

gap> LoadPackage("WPrecog");
gap> K := SymmetricGroup(10);;
gap> H := DihedralGroup(IsPermGroup, 10);;
gap> G := WreathProductProductAction(K, H);;
gap> ri := RecogNode(G);;
gap> G := ri!.Grp;
gap> gens := GeneratorsOfGroup(G);;
gap> L := NrMovedPoints(G);;
gap> delta := 1/10;;
gap> l := 3;;
gap> for i in [1 .. l] do
>       gens := WPR_SimpleSingleComponentBaseStep(gens, L, delta);
>    od; time;
gap> S := Group(gens);;
gap> StructureDescription(S);
gap> Socle(G);; time;

gap> simpleGroupFamily := "Alt";;
gap> origGensS := gens;
gap> hintsForT := rec(upperDegreeBound := 10);
gap> eps := 1/10;
gap> recogData := WPR_RecogniseAlmostSimple(origGensS, simpleGroupFamily, hintsForT, eps);; time;
gap> T := recogData.T;
gap> AutT := recogData.AutT;
gap> slpFuncForT := recogData.slpFuncForT;
gap> stdGensS := recogData.stdGensS;
gap> lambda := recogData.lambda;

gap> domainData := WPR_TopGroupDomain(ri, stdGensS);; time;
gap> t := domainData.t;
gap> domain := domainData.domain;