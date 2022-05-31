 ***,CO
 memory, 250,m
 basis=6-31G
 cartesian
 R1=4.0
 geometry={C;O,C,R1}
 hf
 put,molden,COstretch_631G.molden;
 {FCI;CORE;DUMP}
