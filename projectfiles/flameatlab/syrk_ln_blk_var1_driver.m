% Set the path
% Creating character vector containing a search path that includes all the 
% folders and subfolders in libflameatlab
pathstoadd = genpath('../../libflameatlab');

addpath(pathstoadd);

%% 
% Setting the matrix dimension

m = 100;             
n = 25;
%% 
nb_alg = 10;
% Setting up a triangular matrix
temp= randi( [-3,3], [n,n] );  % random m x m matrix
C = tril(temp) + tril(temp, -1);
 
% Create a random matrix A
A = randi( [-3,3], [n,m] );
%% 
% Check whether syrk_ln_blk_var1( A, C, nb_alg ) computes same as SRYK

if (isequal( syrk_ln_blk_var1( A, C, nb_alg), tril(A*A.' + C)))
    disp( 'All seems well' );
else
    disp( 'Trouble in paradise' )
end
