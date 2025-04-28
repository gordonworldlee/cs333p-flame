
% Copyright 2025 The University of Texas at Austin
%
% For licensing information see
%                http://www.cs.utexas.edu/users/flame/license.html 
%                                                                                 
% Programmed by: gord
%                aniketg@utexas.edu

function [ B_out ] = Trmm_runn_blk_var1( U, B, nb_alg )

  [ UTL, UTR, ...
    UBL, UBR ] = FLA_Part_2x2( U, ...
                               0, 0, 'FLA_BR' );

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_RIGHT' );

  while ( size( UBR, 1 ) < size( U, 1 ) )

    b = min( size( UTL, 1 ), nb_alg );

    [ U00, U01, U02, ...
      U10, U11, U12, ...
      U20, U21, U22 ] = FLA_Repart_2x2_to_3x3( UTL, UTR, ...
                                               UBL, UBR, ...
                                               b, b, 'FLA_TL' );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_LEFT' );

    %------------------------------------------------------------%

    B1 = B1 * U11;
    B1 = B1 + B0 * U01;

    %------------------------------------------------------------%

    [ UTL, UTR, ...
      UBL, UBR ] = FLA_Cont_with_3x3_to_2x2( U00, U01, U02, ...
                                             U10, U11, U12, ...
                                             U20, U21, U22, ...
                                             'FLA_BR' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_RIGHT' );

  end

  B_out = [ BL, BR ];

end
