LoadPackage("WPE", false);
LoadPackage("WPrecog", false);

# Permutation Group : Product Action
Read("examples/difficultPermPrimitiveLargeBase.gi");;

# Recognition
output := RecogniseWreathProduct(RecogNode(G), rec(
    action := "product action",
));;
res := output.res;

timer := Runtime();;
riG := RecogniseGroup(G);;
timer := Runtime() - timer;;
Print("Total Time: ", FormatFloat(timer / 1000.0), " seconds");