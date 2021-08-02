using Base: String, Integer, print_matrix
using DataStructures 


function create_input(file)
    local  x_pos, y_pos, z_pos, natoms, units,atomname
    local method
    local symm
    local occ,closed
    local com,basisdef
    local nel,sym1,mult,nstates,statesx
    local at,ga,ci,typ,type,atoms
    local c,MO,syms,occs,MOnum,MOenergies

    local cc
    local CIvec,confs,c


   statesx=zeros(Int64,2)
    open( "inputs/abinitio.dat") do f

       

        while ! eof(f)

            s=readline(f)
    
            
           
            if contains(s,"GEOMETRY") 
                 s=readline(f)
                 units=s

                 s=readline(f)

                 natoms=parse(Int64,s)

                 println(natoms)

                 x_pos=zeros(Float64,natoms)
                 y_pos=zeros(Float64,natoms)
                 z_pos=zeros(Float64,natoms)
                 atomname=String[]

           
            
                 for i=1:natoms
                     s=readline(f)
                     s=split(s,' ')
                     push!(atomname,s[1])

                     x_pos[i]=parse(Float64,s[2])
                     y_pos[i]=parse(Float64,s[3])
                     z_pos[i]=parse(Float64,s[4])
                 end

            elseif contains(s,"METHOD")
                s=readline(f)
                method=s
               

            elseif contains(s,"STATESYM")
            
                s=readline(f)
               
                sym1=parse(Int64,s)
                        
            elseif contains(s,"SYM")
                s=readline(f)
                if isequal(s,"ON")
                    s=readline(f)
                    symm=s
                else
                    symm="nosym"
                end
            
            elseif contains(s,"BASIS")
                s=readline(f)
                
                basisdef=s

            elseif contains(s,"OCC")
                s=readline(f)
                occ=parse(Int64,s)
               

            elseif contains(s,"CLOSED")
                s=readline(f)
                closed=parse(Int64,s)

            


            elseif contains(s,"COM")
                s=readline(f)
                if isequal(s,"True")
                    com="mass"
                else
                    com="noorient"
                end


            
            elseif contains(s,"NEL")
                s=readline(f)
                nel=parse(Int64,s)
            
           
            
            elseif contains(s,"MULTIPLICITY")
                s=readline(f)
                if isequal(s,"Singlet")
                    mult=0
                elseif isequal("Doublet")
                    mult=1
                elseif isequal("Triplet")
                    mult=2
                end
            
            elseif contains(s,"NSTATES")
                s=readline(f)
                nstates=s
            
            elseif contains(s,"XSCATSTATES")
                s=readline(f)
                s1=split(strip(s))
                statesx[1]=parse(Int64,s1[1])
                statesx[2]=parse(Int64,s1[2])
                
            end
        end

        println(nstates)
    end





    f=open("inputs/molpro.inp","w")

    println(f,"***, scattering calculation in PYXCAT")
    println(f,"gprint,civector,angles=-1,distance=-1
gthresh,twoint=1.0d-13
gthresh,energy=1.0d-7,gradient=1.0d-4
gthresh,thrpun=0.0001
punch,molpro.pun,new
basis=$basisdef
symmetry,$symm;
orient,$com;
$units;
geomtype=xyz;
geometry={
$natoms
")

    for i=1:natoms
        println(f,"$(atomname[i]) $(x_pos[i]) $(y_pos[i]) $(z_pos[i])")
    end 
println(f,"}
hf
")

if isequal(method,"CAS") | isequal(method,"CASSCF")

    println(f,"{multi,failsafe;
maxiter,40;
occ,$occ
closed,$closed
wf,$nel,$sym1,$mult
state,$nstates
pspace,10.0        
orbital,2140.3;
ORBITAL,IGNORE_ERROR;
ciguess,2501.2 
save,ci=2501.2}")



end

println(f,"put, molden, molpro.mld")
close(f)

run(`pwd`)
run(`rm molpro.pun inputs/molpro.out molpro.mld`)

run(`E:/Molpro/bin/molpro.exe -d inputs/ -s inputs/molpro.inp`)

while !isfile("molpro.pun")
    print("waiting")
end 

c=0
realnum=Int64[]
typ=String[]
ga=Float64[]
ci=Float64[]
atoms=Int64[]
MOnum=Int64[]
MO=Float64[]
syms=Float64[]
occs=Float64[]
MOenergies=Float64[]
lang=Float64[]
mang=Float64[]
nang=Float64[]
cc=0

open( "molpro.mld") do f
    s=readline(f)
    while ! eof(f)

        s=readline(f)

        if contains(s,"[GTO]")

            s=readline(f)
            while ! contains(s,"[MO]")
                if length(strip(s))!=0
                    s1=split(strip(s)," ")
                    s1=[s for s in s1 if !isempty(s)]
                    
                    if ! occursin(r"[a-z]+",s1[1]) && ! occursin(r"[.]+",s1[1])
                        at=s1[1]
                      
                    elseif occursin(r"[a-z]+",s1[1])
                        type=s1[1]
                        ncont=parse(Int64,s1[2])
                  

                        if type=="s"
                            
                         
                            for ii=1:ncont 
                                push!(realnum,cc+1)
                                s=readline(f)
                                s1=split(strip(s)," ") 
                                s1=[s for s in s1 if !isempty(s)]                       
                                push!(typ,type)
                                push!(ga,parse(Float64,replace(s1[1],"D" => "E")))
                                push!(ci,parse(Float64,replace(s1[2],"D" => "E")))
                                push!(atoms,parse(Int64,at))
                                push!(lang,0); push!(mang,0); push!(nang,0);
                            end
                            cc+=1
                        elseif type=="p"
                            
                            
                            for ii=1:ncont
                                s=readline(f)
                                s1=split(strip(s)," ") 
                                s1=[s for s in s1 if !isempty(s)] 
                                push!(realnum,cc+1)      
                                push!(ga,parse(Float64,replace(s1[1],"D" => "E")))
                                push!(ci,parse(Float64,replace(s1[2],"D" => "E")))
                                push!(atoms,parse(Int64,at))
                                push!(lang,1); push!(mang,0); push!(nang,0);
                                push!(realnum,cc+2)
                                push!(ga,parse(Float64,replace(s1[1],"D" => "E")))
                                push!(ci,parse(Float64,replace(s1[2],"D" => "E")))
                                push!(atoms,parse(Int64,at))
                                push!(lang,0); push!(mang,1); push!(nang,0);
                                push!(realnum,cc+3)
                                push!(ga,parse(Float64,replace(s1[1],"D" => "E")))
                                push!(ci,parse(Float64,replace(s1[2],"D" => "E")))
                                push!(atoms,parse(Int64,at))
                                push!(lang,0); push!(mang,0); push!(nang,1);
                            end
                            cc+=3
                        end

                        
                        
                    end
                end
                s=readline(f)            
            end
        end
    
        if contains(s, "[MO]")
            
            s=readline(f)
            
            while ! eof(f) & ! contains(s, "[")
                s=readline(f)
                if contains(s,"Sym") 
                   
                    push!(syms,parse(Float64,match(r"\d*\.?\d*$",s,1).match))
               
                elseif contains(s,"Occ")
                            
                    push!(occs,parse(Float64,match(r"\d*\.?\d*$",s,1).match))
                elseif contains(s,"Ene")
                
                            
                    push!(MOenergies,parse(Float64,match(r"\d*\.?\d*$",s,1).match)) 
                elseif ! contains(s, "Spin") & ! isempty(s)
                   
                               
                    push!(MO,parse(Float64,match(r"\d*\.?\d*$",s,1).match))
                   
                    push!(MOnum,parse(Int64,match(r"(?<!\.)\b[0-9]+\b(?!\.)",s).match)) 
                end

            end
        end        
    end
end
close(f)



maxMO= maximum(MOnum)



countsMO=counter(MOnum)



basnum=countsMO[maxMO]
#println(basnum)

MO=reshape(MO,(maxMO,basnum))


#Base.print_matrix(IOContext(stdout, :limit => true), MO)
#println()





MO_rec=MO[:,realnum]

#println(size(MO_rec[1,:]), size(ga))
            

#Reorder the M coefficients with symmetry

symidx = sortperm(syms)

MO_rec=MO_rec[symidx,:]

#println(symidx)

#Now construct the density matrix for the states considered, if I==J => Elastic, if I/=J => Inelastic
confs=[]
CIvec=zeros(0)
c=1
open("molpro.pun") do f 
    s=readline(f)
    while ! eof(f) & ! contains(s,"---")
        s=readline(f)
        if startswith(s," ")
            s1=split(strip(s)," ")
            s1=[s for s in s1 if !isempty(s)]

           
            if contains(symm,"nosym")
                push!(confs,s1[1])
                println("uptohere")
            else 
                println("program the bloody function for the symmetries")
            end
            
            if c ==1
                CIvec=[parse(Float64,i) for i in s1[2:end]]'
            else
                CIvec=[CIvec ;[parse(Float64,i) for i in s1[2:end]]']
            end
            println(CIvec)
            c+=1
            
        end 
    end

end
print(CIvec)
print(lang)

end



