% function to collapse and then remove one variable from multi-variable
% data array

% data :    multi-dimensional matrix (eg. [V1 V2 V3]) where each dimension
% represents one variable

% vars2collapse :   [V1 V2 ...], where V1 & V2 are variables to be collapsed


function data = collapse (data, vars2collapse)

    % orig size
    origSize = ndims(data);

    % collapse variables
    for v = vars2collapse
        data = nanmean(data,v);
    end
    
    % permute to remove singletons
    vars = 1:origSize;
    permuteOrder = find(~ismember(vars,vars2collapse));
    for v = vars2collapse
        if v~=origSize % permute won't work if add last variable
            permuteOrder = [permuteOrder v];
        end
    end
    data = permute(data,permuteOrder);
    
end




% EDITED 27 June to maintain non-collapsed singleton dims
% function data = collapse (data, vars2collapse)
% 
%     for v = 1:length(vars2collapse)
%         data = nanmean(data,vars2collapse(v));
%     end
%     data = squeeze(data);
%     
% end
