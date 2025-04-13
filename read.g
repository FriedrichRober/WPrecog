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

ReadPackage( "WPrecog", "gap/Quick.gi");

ReadPackage( "WPrecog", "gap/subroutines/Init.gi");

ReadPackage( "WPrecog", "gap/subroutines/SingleComponentGroup/Helpers.gi");
ReadPackage( "WPrecog", "gap/subroutines/SingleComponentGroup.gi");

ReadPackage( "WPrecog", "gap/subroutines/RecogniseBaseComponent.gi");

ReadPackage( "WPrecog", "gap/subroutines/RecogniseTopComponent.gi");

ReadPackage( "WPrecog", "gap/subroutines/TopGroupDomainAndAction/ViaLinearAction.gi");
ReadPackage( "WPrecog", "gap/subroutines/TopGroupDomainAndAction/ViaProductAction.gi");
ReadPackage( "WPrecog", "gap/subroutines/TopGroupDomainAndAction/ViaSubgroup.gi");
ReadPackage( "WPrecog", "gap/subroutines/TopGroupDomainAndAction.gi");

ReadPackage( "WPrecog", "gap/subroutines/Image/ViaConjugationAction.gi");
ReadPackage( "WPrecog", "gap/subroutines/Image/ViaLinearAction.gi");
ReadPackage( "WPrecog", "gap/subroutines/Image/ViaProductAction.gi");
ReadPackage( "WPrecog", "gap/subroutines/Image.gi");

ReadPackage( "WPrecog", "gap/Recog.gi");
