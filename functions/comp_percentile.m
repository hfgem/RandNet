function x = comp_percentile(dataset,value)
    perc = prctile(dataset,0.01:0.01:100);
    [~,index] = min(abs(perc'-value));
    x = index*0.01;
end