
/* Copyright 2025 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: gord
                  aniketg@utexas.edu
*/

#include "FLAME.h"

int syrk_ln_unb_var1( FLA_Obj A, FLA_Obj C )
{
    FLA_Obj AT,              A0,
            AB,              a1t,
                            A2;

    FLA_Obj CTL,   CTR,      C00,  c01,     C02, 
            CBL,   CBR,      c10t, gamma11, c12t,
                            C20,  c21,     C22;

    FLA_Part_2x1( A,    &AT, 
                        &AB,            0, FLA_BOTTOM );

    FLA_Part_2x2( C,    &CTL, &CTR,
                        &CBL, &CBR,     0, 0, FLA_BR );

    while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

        FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                                &a1t, 
                            /* ** */            /* *** */
                            AB,                &A2,        1, FLA_TOP );

        FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00,  &c01,     /**/ &C02,
                                                    &c10t, &gamma11, /**/ &c12t,
                            /* ************* */   /* ************************** */
                            CBL, /**/ CBR,       &C20,  &c21,     /**/ &C22,
                            1, 1, FLA_TL );

        /*------------------------------------------------------------*/


        // c21 = c21 + A2 * a1t.'
        FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A2, a1t, FLA_ONE, c21 );

        // gamma11 = gamma11 + a1t * a1t' 
        FLA_Dots( FLA_ONE, a1t, a1t, FLA_ONE, gamma11 );
    
        /*------------------------------------------------------------*/

        FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                /* ** */           /* *** */
                                                    a1t, 
                                &AB,                A2,     FLA_BOTTOM );

        FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00,  /**/ c01,     C02,
                                /* ************** */  /* ************************ */
                                                        c10t, /**/ gamma11, c12t,
                                &CBL, /**/ &CBR,       C20,  /**/ c21,     C22,
                                FLA_BR );

    }

    return FLA_SUCCESS;
}

