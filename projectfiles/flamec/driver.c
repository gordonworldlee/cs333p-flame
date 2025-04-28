#include <stdio.h>
#include <math.h>
#include <time.h>
#include "FLAME.h"

/* Function prototypes */
int trmm_runn_unb_var1(FLA_Obj U, FLA_Obj B);
int trmm_runn_blk_var1(FLA_Obj U, FLA_Obj B, int nb_alg);
int syrk_ln_unb_var3( FLA_Obj A, FLA_Obj C );
int syrk_ln_blk_var3( FLA_Obj A, FLA_Obj C, int nb_alg );

int main(int argc, char *argv[])
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
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Uobj);
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Bobj);
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Bold);
        FLA_Obj_create(FLA_DOUBLE, n, n, 1, n, &Bref);

        /* Generate matrices */
        FLA_Random_tri_matrix(FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, Uobj);
        FLA_Random_matrix(Bold);

        /* Time reference implementation */
        dtime_best_ref = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Bold, Bref);
            dtime = FLA_Clock();
            FLA_Trmm(FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                    FLA_NONUNIT_DIAG, FLA_ONE, Uobj, Bref);
            dtime = FLA_Clock() - dtime;
            dtime_best_ref = (irep == 0) ? dtime : fmin(dtime, dtime_best_ref);
        }
        printf("data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best_ref);
        fflush(stdout);

        /* Time unblocked variant */
        dtime_best_unb = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Bold, Bobj);
            dtime = FLA_Clock();
            trmm_runn_unb_var1(Uobj, Bobj);
            dtime = FLA_Clock() - dtime;
            dtime_best_unb = (irep == 0) ? dtime : fmin(dtime, dtime_best_unb);
        }
        diff_unb = FLA_Max_elemwise_diff(Bobj, Bref);
        printf("data_unb_var1( %d, 1:3 ) = [ %d %le %le ];\n",
               i, n, dtime_best_unb, diff_unb);
        fflush(stdout);

        /* Time blocked variant */
        dtime_best_blk = 0.0;
        for (irep = 0; irep < nrepeats; irep++) {
            FLA_Copy(Bold, Bobj);
            dtime = FLA_Clock();
            trmm_runn_blk_var1(Uobj, Bobj, nb_alg);
            dtime = FLA_Clock() - dtime;
            dtime_best_blk = (irep == 0) ? dtime : fmin(dtime, dtime_best_blk);
        }
        diff_blk = FLA_Max_elemwise_diff(Bobj, Bref);
        printf("data_blk_var1( %d, 1:3 ) = [ %d %le %le ];\n",
            i, n, dtime_best_blk, diff_blk);
        fflush(stdout);

        /* Free matrices */
        FLA_Obj_free(&Uobj);
        FLA_Obj_free(&Bobj);
        FLA_Obj_free(&Bref);
        FLA_Obj_free(&Bold);

        i++;
    }

    FLA_Finalize();

    FLA_Obj Cobj, A, Cold, Cref;
    FLA_Init();
    int m = 100;
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
            syrk_ln_unb_var3( A, Cobj );
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
            syrk_ln_blk_var3( A, Cobj, nb_alg );
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
