(******************************************************************)
(*     Restriction file for DipoleDM_full.fr                                                     *)
(*                                                                                                *)                                            
(*     chi1 couples only to the the Dark Scalar and Higgs    *)
(******************************************************************)

M$Restrictions = {
        Caxx0 ->0,
	  	Caxx1->0,
		Caxx2->0,
		Caxx10->0,
		Caxx20->0,
		Caxx21->0,
		Chxx2->0,
		Chxx20->0,
		Chxx21->0,
	    Chxx1->0,
        Chxx0->0,
	    ychi0->0,
		ychi1->0,
		ychi2->0,		
	    ychi10->0,
		Chxx[i_?NumericQ, j_?NumericQ] :> 0 /; (i+j =!= 3),
		Caxx[i_?NumericQ, j_?NumericQ] :> 0,		
		ychi[1,1] :> 0,
		ychi[2,2] :> 0,
		ychi[3,3] :> 0,
		ychi[2,1] :> 0,
		ychi[1,2] :> 0
}
