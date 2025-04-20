function [value] = RMSE_single(x1, GT1, match_range)
    value = (sum((x1(match_range)-GT1(match_range)).^2)/length(match_range))^0.5;
end