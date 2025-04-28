% Set the path
% Creating character vector containing a search path that includes all the 
% folders and subfolders in libflameatlab
pathstoadd = genpath('../../libflameatlab');

addpath(pathstoadd);

%% 
% Setting the matrix dimension

m = 47;              % problem sizes
n = 35;
%% 

% Setting up a triangular matrix
U = randi( [1,3], [n,n] );  % random m x m matrix
U = triu( U );              % make the matrix upper triangular
 
% Create a random matrix B

B = randi( [-3,3], [m,n] );
%% 
% Check whether trmm_runn_unb_var1( B, U ) computes the same as B * U

if ( isequal( trmm_runn_unb_var1( B, U ), B*U ) )
    disp( 'All seems well' );
else
    disp( 'Trouble in paradise' )
end
