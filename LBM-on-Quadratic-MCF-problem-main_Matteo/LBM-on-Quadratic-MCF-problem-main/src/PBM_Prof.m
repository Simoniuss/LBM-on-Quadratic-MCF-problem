function [ x , status ] =  PBM( f , varargin )

% function [ x , status ] = PBM( f , x , mu , m1 , eps , MaxIter , MInf )
%
% Apply the Proximal Bundle Method for the minimization of the provided
% function f, which must have the following interface:
%
%   [ v , g ] = f( x )
%
% Input:
%
% - x is either a [ n x 1 ] real (column) vector denoting the input of
%   f(), or [] (empty).
%
% Output:
%
% - v (real, scalar): if x == [] this is the best known lower bound on
%   the unconstrained global optimum of f(); it can be -Inf if either f()
%   is not bounded below, or no such information is available. If x ~= []
%   then v = f(x).
%
% - g (real, [ n x 1 ] real vector): this also depends on x. if x == []
%   this is the standard starting point from which the algorithm should
%   start, otherwise it is a subgradient of f() at x (possibly the
%   gradient, but you should not apply this algorithm to a differentiable
%   f)
%
% The other [optional] input parameters are:
%
% - x (either [ n x 1 ] real vector or [], default []): starting point.
%   If x == [], the default starting point provided by f() is used.
%
% - mu (real scalar, optional, default value 1): the fixed weight to be
%   given to the stabilizing term throughout all the algorithm. It must
%   be a strictly positive number.
%
% - m1 (real scalar, optional, default value 0.01): parameter of the
%   Armijo-like condition used to declare a Serious Step; has to be in
%   [0,1).
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   direction (optimal solution of the master problem) is less than or
%   equal to mu * eps. If a negative value is provided, this is used in a
%   *relative* stopping criterion: the algorithm is stopped when the norm
%   of the direction is less than or equal to
%      mu * (- eps) * || norm of the first gradient ||.
%
% - MaxIter (integer scalar, optional, default value 100): the maximum
%   number of iterations.
%
% - MInf (real scalar, optional, default value -Inf): if the algorithm
%   determines a value for f() <= MInf this is taken as an indication that
%   the problem is unbounded below and computation is stopped
%   (a "finite -Inf").
%
% Output:
%
% - x ([ n x 1 ] real column vector): the best solution found so far.
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution;
%
%   = 'unbounded': the algorithm has determined an extrenely large negative
%     value for f() that is taken as an indication that the problem is
%     unbounded below (a "finite -Inf", see MInf above)
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the best solution found so far, but not
%     necessarily the optimal one
%
%   = 'error': the solver of the Master Problem (Cplex used through YALMIP,
%     so that it can easily be switched to any other QP solver) reported
%     some error, which requires to stop optimization
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 02-09-21
 Version 0.20
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = 2;
% 1 = the level sets of f and the trajectory are plotted (when n = 2)
% 2 = the function value / gap are plotted
% all the rest: nothing is plotted

if Plotf == 2
   gap = [];
end

Interactive = false;  % if we pause at every iteration

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isa( f , 'function_handle' )
   error( 'f not a function' );
end

if isempty( varargin ) || isempty( varargin{ 1 } )
   [ fStar , x ] = f( [] );
else
   x = varargin{ 1 };
   if ~ isreal( x )
      error( 'x not a real vector' );
   end

   if size( x , 2 ) ~= 1
      error( 'x is not a (column) vector' );
   end
   
   fStar = f( [] );
end

n = size( x , 1 );

if length( varargin ) > 1
   mu = varargin{ 2 };
   if ~ isscalar( mu )
      error( 'mu is not a real scalar' );
   end
   if mu <= 0
      error( 'mu must be > 0' );
   end       
else
   mu = 1;
end

if length( varargin ) > 2
   m1 = varargin{ 3 };
   if ~ isscalar( m1 )
      error( 'm1 is not a real scalar' );
   end
   if m1 < 0 || m1 > 1
      error( 'm1 is not in [ 0 , 1 ]' );
   end       
else
   m1 = 0.01;
end

if length( varargin ) > 3
   eps = varargin{ 4 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 4
   MaxIter = round( varargin{ 5 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 100;
end

if length( varargin ) > 5
   MInf = varargin{ 6 };
   if ~ isscalar( MInf )
      error( 'MInf is not a real scalar' );
   end
else
   MInf = - Inf;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Proximal Bundle method\n');
if fStar > - Inf
   fprintf( 'iter\trel gap\t\tbest gap\t|| d ||\t\tstep\n\n');
else
   fprintf( 'iter\tf(x)\t\tf best\t\t|| d ||\t\tstep\n\n');
end

% compute first function and subgradient- - - - - - - - - - - - - - - - - -

[ fx , g ] = f( x );

fbest = fx;

G = g';           % matrix of subgradients
F = fx - g' * x;  % vector of translated function values
% each ( fxi , gi , xi ) gives the constraint
%
%  v >= fxi + gi' * ( x + d - xi ) = gi' * ( x + d ) + ( fi - gi' * xi )
%
% so we just keep the single constant fi - gi' * xi instead of xi

if eps < 0
   ng0 = - ng;  % norm of first subgradient: why is there a "-"? ;-)
else
   ng0 = 1;     % un-scaled stopping criterion
end

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

iter = 1;
while true

    % construct the master problem- - - - - - - - - - - - - - - - - - - - -
	% PROBLEM 26 OF THE SURVEY
    d = sdpvar( n , 1 );
    v = sdpvar( 1 , 1 );

	%CONSTRAINTS
    M = [ v * ones( size( G , 1 ) , 1 ) >= F + G * ( d + x )  ];

    % this is close to cheating: use information about f_* in the model
    % (but a lower bound on f_* may actually be available in practice)
    %if fStar > - Inf
    %   M = M + [ v >= fStar ];
    %end
    
	%OBJECTIVE
    c = v + mu * norm( d )^2 / 2;
    
    % solve the master problem- - - - - - - - - - - - - - - - - - - - - - -

    % these are for Cplex, but Cplex for MATLAB is buggy ...
    ops = sdpsettings( 'solver' , 'cplex' , 'cplex.threads' , 1 , ...
                       'verbose' , 0 );
    % last resort: MATLAB QP solver
    %ops = sdpsettings( 'solver' , 'QUADPROG' , 'verbose' , 0 );

    diagnostics = optimize( M , c , ops );

    if diagnostics.problem ~= 0
       status = 'error';
       break;
    end

    d = value( d );
    v = value( v );
    nd = norm( d );

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if fStar > - Inf
       gapk = ( fx - fStar ) / max( [ abs( fStar ) 1 ] );
       bstgapk = ( fbest - fStar ) / max( [ abs( fStar ) 1 ] );

       fprintf( '%4d\t%1.4e\t%1.4e\t%1.4e' , iter , gapk , bstgapk , nd );

       if Plotf == 2
          gap( end + 1 ) = gapk;
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          ylim( [ 1e-15 1e+1 ] );
       end
    else
       fprintf( '%4d\t%1.8e\t%1.8e\t\t%1.4e' , iter , fx , fbest , nd );

       if Plotf == 2
          gap( end + 1 ) = v;
          plot( gap , 'Color' , 'k' , 'LineWidth' , 2 );
       end    
    end
    
    if Plotf == 2
       xlim( [ 0 MaxIter ] );
       ax = gca;
       ax.FontSize = 16;
       ax.Position = [ 0.03 0.07 0.95 0.92 ];
       ax.Toolbar.Visible = 'off';
    end

    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if mu * nd <= eps * ng0 %IF MAKE ALMOST NO MOVE WRT TO PREVIUOS ITERATE
       status = 'optimal';
       break;
    end
    
    if iter > MaxIter
       status = 'stopped';
       break;
    end

    % compute function and subgradient- - - - - - - - - - - - - - - - - - -

    [ fd , g ] = f( x + d );

    if fd <= MInf
       status = 'unbounded';
       break;
    end
    
    if fd < fbest
       fbest = fd;
    end

    G = [ G ; g' ];
    F = [ F ; fd - g' * ( x + d ) ];

    % SS / NS decision- - - - - - - - - - - - - - - - - - - - - - - - - - -

    if fd <= fx + m1 * ( v - fx )
       fprintf( '\tSS\n' );

       if n == 2 && Plotf == 1
          PXY = [ x ,  x + d ];
          line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
                'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
                'Color' , [ 0 0 0 ] );
       end
       
       x = x + d;
       fx = fd;
    else
       fprintf( '\tNS\n' );

       if n == 2 && Plotf == 1
          PXY = [ x ,  x + d ];
          line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
                'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
                'Color' , [ 1 0 0 ] );
       end
    end
    
    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    iter = iter + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



