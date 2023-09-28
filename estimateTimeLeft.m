function estimateTimeLeft(ii,jj,N,n)
% ESTIMATETIMELEFT(i,j,N,n)     Estimate the remaining time in a nested
% for-loop, with i looping from 1:N, j looping from 1:n

clc
fprintf(['Currently evaluating: i = ', num2str(ii), ' of ', ...
    num2str(N), '\n'])
if n ~= 1
    fprintf(['   Currently evaluating: j = ', num2str(jj), ' of ', ...
        num2str(n), '\n'])
end
if (ii==1&&jj==1)
    tic
elseif (ii==N&&jj==n)
    clc
else
    telapsed = toc;
    numleft = (N-ii)*n + n-jj+1;
    numdone = N*n-numleft;
    tPerIteration = telapsed/numdone;
    tleft = tPerIteration*numleft;
    tleftmin = floor(tleft/60);
    tleftsec = floor(tleft - tleftmin*60);
    if floor(tleftmin/60)>0
        tlefthour = floor(tleftmin/60);
        tleftmin = floor(tleftmin - tlefthour*60);
        if floor(tlefthour/24) > 0
            tleftday = floor(tlefthour/24);
            tlefthour = floor(tlefthour-tleftday*24);
            fprintf(['Estimated remaining time: ' ...
                num2str(tleftday) 'd '...
                num2str(tlefthour) 'h '...
                num2str(tleftmin) 'min '...
                num2str(tleftsec) 's' '\n'])
        else
            fprintf(['Estimated remaining time: ' ...
                num2str(tlefthour) ' h '...
                num2str(tleftmin) 'min '...
                num2str(tleftsec) 's' '\n'])
        end
    else
        fprintf(['Estimated remaining time: ' num2str(tleftmin) 'min '...
        num2str(tleftsec) 's' '\n'])
    end
end
