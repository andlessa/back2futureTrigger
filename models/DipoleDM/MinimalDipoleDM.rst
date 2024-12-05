(******************************************************************)
(*     Restriction file for DipoleDM_full.fr                                                     *)
(*                                                                                                *)                                            
(*     Only chi1 has couplings to the Dark Scalar    *)
(******************************************************************)

M$Restrictions = {
        Caxx0 ->0,
	  	Caxx1->0,
	    Chxx1->0,
        Chxx0 -> 0,
	    ychi0->0,
	    ychi10->0,
		Chxx[i_?NumericQ, j_?NumericQ] :> 0 /; (i == j),
		Caxx[i_?NumericQ, j_?NumericQ] :> 0 /; (i == j),
		ychi[i_?NumericQ, j_?NumericQ] :> 0 /; (i+j =!= 4)
}
