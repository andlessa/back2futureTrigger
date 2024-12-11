(******************************************************************)
(*     Restriction file for DipoleDM_full.fr                                                     *)
(*                                                                                                *)                                            
(*     chi1 couples only to the the Dark Scalar and Photon    *)
(******************************************************************)

M$Restrictions = {
        Caxx0 ->0,
	  	Caxx1->0,
	    Chxx1->0,
        Chxx0 -> 0,
		Chxx10 -> 0,
	    ychi0->0,
	    ychi10->0,
		Chxx[i_?NumericQ, j_?NumericQ] :> 0,
		Caxx[i_?NumericQ, j_?NumericQ] :> 0 /; (i == j),
		ychi[i_?NumericQ, j_?NumericQ] :> 0 /; (i+j =!= 4)
}
