function [value1,value2,value3,value] = RMSE(x1, GT1, x2, GT2, x3, GT3, match_range)
    value1 = (sum((x1(match_range)-GT1(match_range)).^2)/length(match_range))^0.5;
    value2 = (sum((x2(match_range)-GT2(match_range)).^2)/length(match_range))^0.5;
    value3 = (sum((x3(match_range)-GT3(match_range)).^2)/length(match_range))^0.5;
    value = (value1+value2+value3)/3;
end