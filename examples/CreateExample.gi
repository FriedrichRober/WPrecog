
outputStream := OutputTextFile("examples/difficultExample.gi", false);
SetPrintFormattingStatus(outputStream, false);
AppendTo(outputStream, "G := ", String(G), ";;");
CloseStream(outputStream);

outputStream := OutputTextFile("examples/difficultExample.mag", false);
SetPrintFormattingStatus(outputStream, false);
AppendTo(outputStream, "G := ", ConvertToMagmaInputString(G));
CloseStream(outputStream);