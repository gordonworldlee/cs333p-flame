% Set the path
% Creating character vector containing a search path that includes all the 
% folders and subfolders in libflameatlab
pathstoadd = genpath('../../libflameatlab');

addpath(pathstoadd);

%% 
% Setting the matrix dimension

m = 100;              % problem sizes
n = 25;
nb_alg = 20;
%% 

% Setting up a triangular matrix
U = randi( [1,3], [n,n] );  % random m x m matrix
U = triu( U );              % make the matrix lower triangular
 
% Create a random matrix B

B = randi( [-3,3], [m,n] );
%% 
% Check whether trmm_runn_blk_var1( U, B ) computes the same as L * B

if ( isequal( trmm_runn_blk_var1( U, B, nb_alg ), B*U) )
    disp( 'All seems well' );
else
    disp( 'Trouble in paradise' )
end
