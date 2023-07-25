***, CO
memory, 250,m
basis=6-31G
cartesian
R1=1.50D0
geometry={C;O,C,R1}
hf
put,molden,CO_631G.molden;
{FCI;CORE;DUMP}
