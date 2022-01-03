function x = comp_percentile(dataset,value)
    %______
    %ABOUT: This function takes a dataset and a value as input to calculate
    %the percentile of the value in relation to the dataset.
    %
    %INPUTS:
    %       dataset = data vector or array of values
    %       value = single value
    %
    %OUTPUTS:
    %       x = percentile value rounded to the 2nd decimal place
    %______
    perc = prctile(dataset,0.01:0.01:100);
    [~,index] = min(abs(perc'-value));
    x = index*0.01;
end