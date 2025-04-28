
/* Copyright 2025 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: gord
                  aniketg@utexas.edu
*/

#include "FLAME.h"

int syrk_ln_blk_var3( FLA_Obj A, FLA_Obj C, int nb_alg )
{
    FLA_Obj AT,              A0,
            AB,              A1,
                            A2;

    FLA_Obj CTL,   CTR,      C00, C01, C02, 
            CBL,   CBR,      C10, C11, C12,
                            C20, C21, C22;

    int b;

    FLA_Part_2x1( A,    &AT, 
                        &AB,            0, FLA_BOTTOM );

    FLA_Part_2x2( C,    &CTL, &CTR,
                        &CBL, &CBR,     0, 0, FLA_BR );

    while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

        b = min( FLA_Obj_length( AT ), nb_alg );

        FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                                &A1, 
                            /* ** */            /* ** */
                            AB,                &A2,        b, FLA_TOP );

        FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, &C01, /**/ &C02,
                                                    &C10, &C11, /**/ &C12,
                            /* ************* */   /* ******************** */
                            CBL, /**/ CBR,       &C20, &C21, /**/ &C22,
                            b, b, FLA_TL );

        /*------------------------------------------------------------*/

        FLA_Syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
            FLA_ONE, A1, FLA_ONE, C11 );

        FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A2, A1, FLA_ONE, C21 );

        /*------------------------------------------------------------*/

        FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                /* ** */           /* ** */
                                                    A1, 
                                &AB,                A2,     FLA_BOTTOM );

        FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, /**/ C01, C02,
                                /* ************** */  /* ****************** */
                                                        C10, /**/ C11, C12,
                                &CBL, /**/ &CBR,       C20, /**/ C21, C22,
                                FLA_BR );

    }

    return FLA_SUCCESS;
}

