#
# WPrecog: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# This file is a script which compiles the package manual.
#
if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;

AutoDoc( rec( scaffold := rec(
        # bib := "wpr",
        includes := [
            "intro.xml",
            "functions.xml",
            ],
        ),
        extract_examples := true,
        autodoc := false ) );

Exec("dev/tests_doc/processTests.sh");