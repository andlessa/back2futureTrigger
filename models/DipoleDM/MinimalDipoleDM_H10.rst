(******************************************************************)
(*     Restriction file for DipoleDM_full.fr                                                     *)
(*                                                                                                *)                                            
(*     chi1 couples only to the the Dark Scalar and Higgs    *)
(******************************************************************)

M$Restrictions = {
        Caxx0 ->0,
	  	Caxx1->0,
		Caxx10->0,
	    Chxx1->0,
        Chxx0 -> 0,
	    ychi0->0,
		Chxx[i_?NumericQ, j_?NumericQ] :> 0 /; (i == j),
		Caxx[i_?NumericQ, j_?NumericQ] :> 0,
		ychi[i_?NumericQ, j_?NumericQ] :> 0 /; (i+j < 3)
}
