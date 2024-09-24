function y = DS1054Z_Convert(data_obj)
    y = (single(data_obj.rawData)-data_obj.yReference)*data_obj.yIncrement;
end