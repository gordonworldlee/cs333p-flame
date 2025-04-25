#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void trmm_runn_unb_var1( FLA_Obj, FLA_Obj );
int main(int argc, char *argv[])
{
    int n, nfirst, nlast, ninc, i, irep, nrepeats;

    double
        dtime, dtime_best, 
        diff;

    dtime_best = 0.0;

    FLA_Obj
        Uobj, Bobj, Bold, Bref;

    /* Initialize FLAME. */
    FLA_Init( );

    /* Every time trial is repeated "repeat" times and the fastest run in recorded */
    printf( "%% number of repeats:" );
    scanf( "%d", &nrepeats );
    printf( "%% %d\n", nrepeats );

    /* Timing trials for matrix sizes n=nfirst to nlast in increments 
       of ninc will be performed. */
    printf( "%% enter nfirst, nlast, ninc:" );
    scanf( "%d%d%d", &nfirst, &nlast, &ninc );
    printf( "%% %d %d %d \n", nfirst, nlast, ninc );
    fflush( stdout );

    i = 1;
    for ( n=nfirst; n<= nlast; n+=ninc ){

        /* Allocate space for the matrices and vectors */
        FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Uobj );
        FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bobj );
        FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bold );
        FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bref );

        /* Generate random matrix L and B */
        //FLA_Random_matrix( Uobj );
        FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, Uobj );
        FLA_Random_matrix( Bold );



        for ( irep=0; irep<nrepeats; irep++ ) {
            /* Time reference implementation (from libflame) */
            FLA_Copy( Bold, Bref );

            /* start clock */
            dtime = FLA_Clock();

            /* Compute Bref = Bref * U where U is upper triangular... */
            FLA_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, Uobj, Bref );

            /* stop clock */
            dtime = FLA_Clock() - dtime;

            if ( irep == 0 ) 
                dtime_best = dtime;
            else
                dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
        }

        printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
        fflush( stdout );

        /* Time your unblocked Variant 1 */

        for ( irep=0; irep<nrepeats; irep++ ){
            /* Copy vector yold to y */
            FLA_Copy( Bold, Bobj );

            /* start clock */
            dtime = FLA_Clock();

            /* Comment out the below call and call your routine instead */
            trmm_runn_unb_var1( Uobj, Bobj );

            /* stop clock */
            dtime = FLA_Clock() - dtime;

            if ( irep == 0 ) 
                dtime_best = dtime;
            else
                dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
        }

        diff = FLA_Max_elemwise_diff( Bobj, Bref );

        printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
                dtime_best, diff  );

        fflush( stdout );

        FLA_Obj_free( &Uobj );
        FLA_Obj_free( &Bobj );
        FLA_Obj_free( &Bref );
        FLA_Obj_free( &Bold );

        i++;
    }
    FLA_Finalize( );

    exit( 0 );
}
