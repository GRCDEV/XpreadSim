function trace2 = ONETrace_extract_from_interval(trace,t1,t2)

function trace2 = remove_duplicate_times(trace) 
    trace = sortrows(trace,2);
    trace(:,1) = floor(trace(:,1));
    jj = 1;
    trace2=trace;
    trace2(1,:) = trace(1,:);
    for ii = 2:length(trace)
        if trace(ii,1) == trace(ii-1,1) && trace(ii,2) == trace(ii-1,2)
            trace2(jj,:) = trace(ii,:);
        else
            trace2(jj,:) = trace(ii,:);
            jj = jj + 1;
        end

    end
    trace2 = trace2(1:jj-1,:);
    trace2 = sortrows(trace2,1);
end

    % Extract ONE trace from t1 to t2
    for i = 1:length(trace)
        if trace(i,1) >= t1
            break;
        end
    end

    iStart = i;

    for i = iStart:length(trace)
        if trace(i,1) >= t2
            break;
        end
    end
    
    iEnd = i;

    trace2 = trace(iStart:iEnd,:);
        
    trace2 = remove_duplicate_times(trace2);
end
