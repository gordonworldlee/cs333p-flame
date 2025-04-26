#include <stdio.h>
#include <math.h>
#include <time.h>
#include "FLAME.h"

/* Function prototypes */
int syrk_ln_unb_var1( FLA_Obj A, FLA_Obj C );
int Syrk_ln_blk_var1( FLA_Obj A, FLA_Obj C, int nb_alg );

int test_syrk(int argc, char *argv[])
{
    int n, nfirst, nlast, ninc, i, irep, nrepeats;
    double dtime, dtime_best_ref, dtime_best_unb, dtime_best_blk;
    double diff_unb, diff_blk;
    // Block size for blocked algorithm
    const int nb_alg = 100; 

    FLA_Obj Uobj, Bobj, Bold, Bref;

    FLA_Init();

    /* Input parameters */
    printf("%% number of repeats:");
    scanf("%d", &nrepeats);
    printf("%% %d\n", nrepeats);

    printf("%% enter nfirst, nlast, ninc:");
    scanf("%d%d%d", &nfirst, &nlast, &ninc);
    printf("%% %d %d %d \n", nfirst, nlast, ninc);
    fflush(stdout);

    i = 1;
    for (n = nfirst; n <= nlast; n += ninc) {
        /* Matrix allocation */
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cobj);
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cold);
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Cref);
        FLA_Obj_create(FLA_DOUBLE, n, m, 1, n, &A);
        /* Generate matrices */
        FLA_Random_symm_matrix( FLA_LOWER_TRIANGULAR, Cold );

        /* Time reference implementation */
        dtime_best_ref = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Cold, Cref);
            dtime = FLA_Clock();
            FLA_Syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                FLA_ONE, A, FLA_ONE, Cref );
            dtime = FLA_Clock() - dtime;
            dtime_best_ref = (irep == 0) ? dtime : fmin(dtime, dtime_best_ref);
        }
        printf("data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best_ref);
        fflush(stdout);

        /* Time unblocked variant */
        dtime_best_unb = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Cold, Cobj);
            dtime = FLA_Clock();
            Syrk_ln_unb_var1( A, Cobj );
            dtime = FLA_Clock() - dtime;
            dtime_best_unb = (irep == 0) ? dtime : fmin(dtime, dtime_best_unb);
        }
        diff_unb = FLA_Max_elemwise_diff(Cobj, Cref);
        printf("data_unb_var1( %d, 1:3 ) = [ %d %le %le ];\n",
               i, n, dtime_best_unb, diff_unb);
        fflush(stdout);

        /* Time blocked variant */
        dtime_best_blk = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Cold, Cobj);
            dtime = FLA_Clock();
            Syrk_ln_blk_var1( A, Cobj, nb_alg );
            dtime = FLA_Clock() - dtime;
            dtime_best_blk = (irep == 0) ? dtime : fmin(dtime, dtime_best_blk);
        }
        diff_blk = FLA_Max_elemwise_diff(Cobj, Cref);
        printf("data_blk_var1( %d, 1:3 ) = [ %d %le %le ];\n",
            i, n, dtime_best_blk, diff_blk);
        fflush(stdout);

        /* Free matrices */
        FLA_Obj_free(&A);
        FLA_Obj_free(&Cobj);
        FLA_Obj_free(&Cref);
        FLA_Obj_free(&Cold);

        i++;
    }

    FLA_Finalize();
    return 0;
}
