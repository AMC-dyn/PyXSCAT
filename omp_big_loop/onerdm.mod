  `7  _   k820309    ,          2021.5.0    �ٿd                                                                                                          
       one_rdm.f90 ONERDM                                                     
       #         @                                                      #CONFS    #EP2    #NDIFF              
                                                                 &                   &                                                                                                                     &                   &                                                                                                                     &                   &                                           (        D                                                   :                    #I              &                                                     
                                            %         @                                                           #DET1 	   #DET2    #NINT 
            
                                 	                          p        4 5 � p        r 
   p          4 5 � p        r 
     p            4 5 � p        r 
     p                                   
                                                           p        4 5 � p        r 
   p          4 5 � p        r 
     p            4 5 � p        r 
     p                                    
                                 
           #         @                                                       #CONFS    #CIVS    #ONERDM    #MAXNMO    #STATE1    #STATE2              
  @                                                               &                   &                                                     
  @                                                 
              &                   &                                                   D @                                                 
               &                   &                                                     D @                                                    
  @                                                   
  @                                         #         @                                                      #CONFS    #CIVS    #NDIFF    #EP2    #ONERDM    #MAXNMO    #STATE1    #STATE2              
                                                                 &                   &                                                     
                                                    
              &                   &                                                     
                                                     
             &                   &                                                     
                                                     	             &                   &                                                   D                                                  
               &                   &                                                     D                                                      
                                                      
                                            #         @                                                       #FILE_READ    #ONERDM    #MAXNMO                                           
                                                                    D                                                  
               &                   &                                                     
                                            #         @                                                        #FILE_READ !   #NUMBERLINES "   #NEWDAT #   #IREP $   #START %   #END &   #MAT '   #TOTAL (                                                                                      
                                 !                                     
                                 "                     
                                #                                 &                   &                                                     
 @                              $                                 &                                                     
                                 %                                 &                                                     
                                 &                                 &                                                   D                                '                                  &                   &                                                   D                                (                   
               &                                           %         @�                               )                           #N *   #K +             
  @                              *                     
  @                              +           #         @                                   ,                   #CREATETWORDM_BIT%PHASE_DBL -   #FILE_READ .   #LENGTH /   #MAT 0   #TOTAL 1                    �                           -                   
                                                    T
W
p          n
       
                       �?        1.D0  n
         
                        �          h  p          p           & p         p            p                                                            
                                 .                                                                     /                    D                                0                    ;              &                   &                                                   D                                1                   
 :              &                                           #         @                                   2                    #MAT 3   #TWORDM 4   #ONERDM 5   #NEL 6             
 @                              3                    E             &                   &                                                     
                                  4                   
 F             &                                                   D                                 5                   
 G              &                   &                                                     
                                 6           #         @                                   7                    #FILE_READ 8   #ONERDM 9   #MAXNMO :   #NUMBERLINES ;   #NEWDAT <   #IREP =                                                                    
                                 8                                   D                               9                   
 K              &                   &                                                     
D                                :                                                      ;                      
                                 <                    O             &                   &                                                     
                                 =                    H             &                                           #         @                                   >                    #DET ?   #NDET A   #COEF B   #MO_NUM C   #NINT @   #DENSITY_MATRIX D            
  @                              ?                     T       p        p        p        5 � p        r @   p          5 � p        r @     p          5 � p        r A       5 � p        r @     p          5 � p        r A                               
                        @         A                    
                                  B                    
 U   p          5 � p        r A       5 � p        r A                               
                                  C                     
  @                     @         @                    D                                 D                    
 V      p        5 � p        r C   p          5 � p        r C     5 � p        r C       5 � p        r C     5 � p        r C                     #         @                                  E                    #DET1 F   #DET2 H   #EXC I   #DEGREE J   #PHASE K   #NINT G            
  @                              F                     `     p        4 5 � p        r G   p          4 5 � p        r G     p            4 5 � p        r G     p                                   
  @                              H                     a     p        4 5 � p        r G   p          4 5 � p        r G     p            4 5 � p        r G     p                                    D @   �                           I                    b    p 	         & p         p          p          p            p          p          p                                    D                                J                      D @                               K     
                 
  @                              G           #         @                                  L                   #GET_DOUBLE_EXCITATION%PHASE_DBLE M   #DET1 N   #DET2 P   #EXC Q   #PHASE R   #NINT O                    �                            M                   
                                                    T
W
p          n
       
                       �?        1.D0  n
         
                        �          h  p          p           & p         p            p                                                           
                                 N                     X     p        4 5 � p        r O   p          4 5 � p        r O     p            4 5 � p        r O     p                                   
                                 P                     Y     p        4 5 � p        r O   p          4 5 � p        r O     p            4 5 � p        r O     p                                    D     �                           Q                    Z    p 	         & p         p          p          p            p          p          p                                    D                                 R     
                 
                                 O           #         @                                  S                   #GET_SINGLE_EXCITATION%PHASE_DBLE T   #DET1 U   #DET2 W   #EXC X   #PHASE Y   #NINT V                    �                            T                   
                                                    T
W
p          n
       
                       �?        1.D0  n
         
                        �          h  p          p           & p         p            p                                                           
                                 U                     \     p        4 5 � p        r V   p          4 5 � p        r V     p            4 5 � p        r V     p                                   
                                 W                     ]     p        4 5 � p        r V   p          4 5 � p        r V     p            4 5 � p        r V     p                                    D     �                           X                    ^    p 	         & p         p          p          p            p          p          p                                    D                                 Y     
                 
                                 V           #         @                                   Z                   #MCCI_TO_BIT%PHASE_DBL [   #FILE_READ \   #FILE_WRITE ]   #NUMBERLINES ^                    �                            [                   
                                                    T
W
p          n
       
                       �?        1.D0  n
         
                        �          h  p          p           & p         p            p                                                            
                                 \     <                                
                                 ]     <                                
                                 ^              �         fn#fn    �   @   J   TWORDMS '   �   g       MAXCOINCIDENCE+TWORDMS -   b  �   a   MAXCOINCIDENCE%CONFS+TWORDMS +     �   a   MAXCOINCIDENCE%EP2+TWORDMS -   �  �   a   MAXCOINCIDENCE%NDIFF+TWORDMS ,   N  �       INTEGER2BINARY_ORBS+TWORDMS .   �  @   a   INTEGER2BINARY_ORBS%I+TWORDMS    1  n       EXCITATIONS !   �    a   EXCITATIONS%DET1 !   �    a   EXCITATIONS%DET2 !   �  @   a   EXCITATIONS%NINT    �  �       ONERDM_CREAT #   �  �   a   ONERDM_CREAT%CONFS "   0  �   a   ONERDM_CREAT%CIVS $   �  �   a   ONERDM_CREAT%ONERDM $   x	  @   a   ONERDM_CREAT%MAXNMO $   �	  @   a   ONERDM_CREAT%STATE1 $   �	  @   a   ONERDM_CREAT%STATE2    8
  �       CREATEONERDM #   �
  �   a   CREATEONERDM%CONFS "   }  �   a   CREATEONERDM%CIVS #   !  �   a   CREATEONERDM%NDIFF !   �  �   a   CREATEONERDM%EP2 $   i  �   a   CREATEONERDM%ONERDM $     @   a   CREATEONERDM%MAXNMO $   M  @   a   CREATEONERDM%STATE1 $   �  @   a   CREATEONERDM%STATE2    �  �       ONE_RDM_SLOW '   Y  P   a   ONE_RDM_SLOW%FILE_READ $   �  �   a   ONE_RDM_SLOW%ONERDM $   M  @   a   ONE_RDM_SLOW%MAXNMO %   �  �       CREATETWORDM_BIT_FCI /   |  P   a   CREATETWORDM_BIT_FCI%FILE_READ 1   �  @   a   CREATETWORDM_BIT_FCI%NUMBERLINES ,     �   a   CREATETWORDM_BIT_FCI%NEWDAT *   �  �   a   CREATETWORDM_BIT_FCI%IREP +   <  �   a   CREATETWORDM_BIT_FCI%START )   �  �   a   CREATETWORDM_BIT_FCI%END )   T  �   a   CREATETWORDM_BIT_FCI%MAT +   �  �   a   CREATETWORDM_BIT_FCI%TOTAL    �  ^       COMBO    �  @   a   COMBO%N    "  @   a   COMBO%K !   b  �       CREATETWORDM_BIT +   �  �     CREATETWORDM_BIT%PHASE_DBL +   �  P   a   CREATETWORDM_BIT%FILE_READ (   �  @   a   CREATETWORDM_BIT%LENGTH %     �   a   CREATETWORDM_BIT%MAT '   �  �   a   CREATETWORDM_BIT%TOTAL     A  r       ONE_RDM_TWO_RDM $   �  �   a   ONE_RDM_TWO_RDM%MAT '   W  �   a   ONE_RDM_TWO_RDM%TWORDM '   �  �   a   ONE_RDM_TWO_RDM%ONERDM $   �  @   a   ONE_RDM_TWO_RDM%NEL    �  �       ONE_RDM_BIT &   �  P   a   ONE_RDM_BIT%FILE_READ #   �  �   a   ONE_RDM_BIT%ONERDM #   �  @   a   ONE_RDM_BIT%MAXNMO (   �  @   a   ONE_RDM_BIT%NUMBERLINES #     �   a   ONE_RDM_BIT%NEWDAT !   �  �   a   ONE_RDM_BIT%IREP '   8   �       COMPUTE_DENSITY_MATRIX +   �   d  a   COMPUTE_DENSITY_MATRIX%DET ,   +"  @   a   COMPUTE_DENSITY_MATRIX%NDET ,   k"  �   a   COMPUTE_DENSITY_MATRIX%COEF .   #  @   a   COMPUTE_DENSITY_MATRIX%MO_NUM ,   _#  @   a   COMPUTE_DENSITY_MATRIX%NINT 6   �#  $  a   COMPUTE_DENSITY_MATRIX%DENSITY_MATRIX    �$  �       GET_EXCITATION $   I%    a   GET_EXCITATION%DET1 $   Y&    a   GET_EXCITATION%DET2 #   i'  �   a   GET_EXCITATION%EXC &   M(  @   a   GET_EXCITATION%DEGREE %   �(  @   a   GET_EXCITATION%PHASE $   �(  @   a   GET_EXCITATION%NINT &   )  �       GET_DOUBLE_EXCITATION 1   �)  �     GET_DOUBLE_EXCITATION%PHASE_DBLE +   5+    a   GET_DOUBLE_EXCITATION%DET1 +   E,    a   GET_DOUBLE_EXCITATION%DET2 *   U-  �   a   GET_DOUBLE_EXCITATION%EXC ,   9.  @   a   GET_DOUBLE_EXCITATION%PHASE +   y.  @   a   GET_DOUBLE_EXCITATION%NINT &   �.  �       GET_SINGLE_EXCITATION 1   Y/  �     GET_SINGLE_EXCITATION%PHASE_DBLE +   �0    a   GET_SINGLE_EXCITATION%DET1 +   �1    a   GET_SINGLE_EXCITATION%DET2 *   3  �   a   GET_SINGLE_EXCITATION%EXC ,   �3  @   a   GET_SINGLE_EXCITATION%PHASE +   %4  @   a   GET_SINGLE_EXCITATION%NINT    e4  �       MCCI_TO_BIT &   �4  �     MCCI_TO_BIT%PHASE_DBL &   �6  P   a   MCCI_TO_BIT%FILE_READ '   �6  P   a   MCCI_TO_BIT%FILE_WRITE (    7  @   a   MCCI_TO_BIT%NUMBERLINES 