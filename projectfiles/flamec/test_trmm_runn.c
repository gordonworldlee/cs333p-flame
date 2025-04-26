#include <stdio.h>
#include <math.h>
#include <time.h>
#include "FLAME.h"
#include "trmm_tests.h"

/* External prototypes for your TRMM implementations */
extern int trmm_runn_unb_var1(FLA_Obj U, FLA_Obj B);
extern int trmm_runn_blk_var1(FLA_Obj U, FLA_Obj B, int nb_alg);

void run_trmm_tests(int nrepeats, int nfirst, int nlast, int ninc, int nb_alg) {
    int i = 1, n, irep;
    double dtime, dtime_best_ref, dtime_best_unb, dtime_best_blk;
    double diff_unb, diff_blk;
    FLA_Obj Uobj, Bobj, Bold, Bref;

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

        /* Free matrices */
        FLA_Obj_free(&Uobj);
        FLA_Obj_free(&Bobj);
        FLA_Obj_free(&Bref);
        FLA_Obj_free(&Bold);

        i++;
    }
}
