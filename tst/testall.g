#
# WPrecog: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "WPrecog" );

docTests := TestDirectory(DirectoriesPackageLibrary( "WPrecog", "tst/files/doc"),
  rec(
    testOptions := rec(
      width := 120,
      compareFunction := "uptowhitespace",

    ),
  )
);

# machineTests := TestDirectory(DirectoriesPackageLibrary( "WPrecog", "tst/files/machine-generated"),
#   rec()
# );

# humanTests := TestDirectory(DirectoriesPackageLibrary( "WPrecog", "tst/files/human-created"),
#   rec()
# );

# if not (docTests and machineTests and humanTests) then
if not (docTests) then
  FORCE_QUIT_GAP(1); # if we ever get here, there was an error
else
  FORCE_QUIT_GAP(0); # everything is fine
fi;
