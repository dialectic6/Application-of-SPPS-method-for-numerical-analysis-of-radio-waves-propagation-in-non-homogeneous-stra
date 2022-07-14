function result = ninteg( dat, bminusa )

% Usage: ninteg(dat, intervallength ).
% Integrates list dat of n=5k+1 data points. Returns list of n points.
%  Uses 5-point integrated interpolating polynomial,
%    applied to points 1:6, 6:10, 11:16, etc.
%  Supposes points equally spaced on interval.
%
n = length(dat);                           % dat is to be integrated
if n < 6
    error( ['ninteg: requires at least 6 data points, received ' num2str(n) ] )
end
n2 = mod(n-1,5);            % number of points to process at the beginning
n1 = n - n2;                % number of points to process afterwards

intmat = [  [475,1427,-798,482,-173,27]/1440    % Formula for
    [28,129,14,14,-6,1]/90              % numerical
    3*[17,73,38,38,-7,1]/160              % integration
    2*[7,32,12,32,7,0]/45                 % A 6 x 5 matrix
    5*[19,75,50,50,75,19]/288
    ] * bminusa/ (n-1);                    % Divide by interval length

inval = 0;
if n2 > 0
    mdat = reshape( dat(1:5), 5, 1 );
    row_n = mdat(1,:);                         % Create final row, so each column
    row_n = [ row_n(2:length(row_n)), dat(6) ];% ends with start of next column
    mdat = [ mdat', row_n' ]';                % Annex the row
    m1 = intmat * mdat;       % numerical integration, gives 5 values
    incr = m1(5,:);           % prepare to calculate cumulative sums
    incr = cumsum( [ 0, incr(1:length(incr)-1) ] );
    incr = [ incr; incr; incr; incr; incr ];
    m1 = m1 + incr;           % matrix now contains the integrals
    result = reshape(m1,1,5);   % convert to a single list of length n
    tresult = result(1:n2);             % take the first n2 values
    inval = tresult(end);
end
% proceed the first part of the list
%if mod(n,5) == 1                           % Verify n is 5k+1
    mdat = reshape( dat(n2+1:n-1), 5, (n1-1)/5 ); % Break into 6 rows: first 5 here
    row_n = mdat(1,:);                         % Create final row, so each column
    row_n = [ row_n(2:length(row_n)), dat(n) ];% ends with start of next column
    mdat = [ mdat', row_n' ]';                % Annex the row
    m1 = intmat * mdat;       % numerical integration, gives 5 values
    incr = m1(5,:);           % prepare to calculate cumulative sums
    incr = cumsum( [ 0, incr(1:length(incr)-1) ] );
    incr = [ incr; incr; incr; incr; incr ];
    m1 = m1 + incr + inval;           % matrix now contains the integrals
    if n2 > 0 
        result = [0, tresult, reshape(m1,1,n1-1) ];   % convert to a single list of length n
    else
        result = [0, reshape(m1,n2+1,n-1) ];   % convert to a single list of length n
    end
% else
%     error( ['ninteg: requires 5k+1 data points, received ' num2str(n) ] )
% end

end
