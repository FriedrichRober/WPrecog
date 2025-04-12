#
# WPrecog: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# Reading the implementation part of the package.
#

# Level 1: Show Current Step in Main Recognition Method
# Level 2: Show Important Bounds in Submethods
# Level 3: Show Progress of Iterations
BindGlobal("WPR_Info", NewInfoClass("WPR_Info"));
SetInfoLevel(WPR_Info, 3);

ReadPackage( "WPrecog", "gap/Random.gi");
ReadPackage( "WPrecog", "gap/SingleComponentGroup.gi");
ReadPackage( "WPrecog", "gap/TopGroupAction.gi");
ReadPackage( "WPrecog", "gap/RecogSubroutines.gi");
ReadPackage( "WPrecog", "gap/Recog.gi");
