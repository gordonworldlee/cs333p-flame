% Set the path
% Creating character vector containing a search path that includes all the 
% folders and subfolders in libflameatlab
pathstoadd = genpath('../../libflameatlab');

addpath(pathstoadd);

%% 
% Setting the matrix dimension

m = 6;              % problem sizes
%% 

% Setting up a triangular matrix
L = randi( [1,3], [m,m] );  % random m x m matrix
L = tril( L );              % make the matrix lower triangular
 
% Create a random matrix B

B = randi( [-3,3], [m,n] );
%% 
% Check whether trmm_runn_unb_var1( L, B ) computes the same as L * B

if ( isequal( trmm_runn_unb_var1( L, B ), L*B ) )
    disp( 'All seems well' );
else
    disp( 'Trouble in paradise' )
end
