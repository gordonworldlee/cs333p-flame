% Set the path
% Creating character vector containing a search path that includes all the 
% folders and subfolders in libflameatlab
pathstoadd = genpath('../../libflameatlab');

addpath(pathstoadd);

%% 
% Setting the matrix dimension

m = 34;             
n = 47;
%% 

% Setting up a triangular matrix
temp= randi( [-3,3], [n,n] );  % random m x m matrix
C = tril(temp) + tril(temp, -1); % Make symmetric matrix
 
% Create a random matrix A
A = randi( [-3,3], [n,m] );
%% 
% Check whether syrk_ln_unb_var3( A, C ) computes the same as C = AA^T + C

if (isequal( syrk_ln_unb_var3( A, C), tril(A*A.' + C)))
    disp( 'All seems well' );
else
    disp( 'Trouble in paradise' )
end
