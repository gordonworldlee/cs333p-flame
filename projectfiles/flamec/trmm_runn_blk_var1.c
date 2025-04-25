
/* Copyright 2025 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: gord
                  aniketg@utexas.edu
*/

#include "FLAME.h"

int trmm_runn_blk_var1( FLA_Obj U, FLA_Obj B, int nb_alg )
{
    FLA_Obj UTL,   UTR,      U00, U01, U02, 
            UBL,   UBR,      U10, U11, U12,
                            U20, U21, U22;

    FLA_Obj BL,    BR,       B0,  B1,  B2;

    int b;

    FLA_Part_2x2( U,    &UTL, &UTR,
                        &UBL, &UBR,     0, 0, FLA_BR );

    FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

    while ( FLA_Obj_length( UBR ) < FLA_Obj_length( U ) ){

        b = min( FLA_Obj_length( UTL ), nb_alg );

        FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, &U01, /**/ &U02,
                                                    &U10, &U11, /**/ &U12,
                            /* ************* */   /* ******************** */
                            UBL, /**/ UBR,       &U20, &U21, /**/ &U22,
                            b, b, FLA_TL );

        FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &B1, /**/ &B2,
                            b, FLA_LEFT );

        /*------------------------------------------------------------*/

        FLA_Trmm(FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, 
            FLA_NONUNIT_DIAG, FLA_ONE, U11, B1);
        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, B0, U01, FLA_ONE, B1);
        FLA_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, B2, U21, FLA_ONE, B1);

        /*------------------------------------------------------------*/

        FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, /**/ U01, U02,
                                /* ************** */  /* ****************** */
                                                        U10, /**/ U11, U12,
                                &UBL, /**/ &UBR,       U20, /**/ U21, U22,
                                FLA_BR );

        FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ B1, B2,
                                FLA_RIGHT );

    }

    return FLA_SUCCESS;
}

