maxDemand_vec = [4000,6000,8000,10000];
% maxDemand_vec = [3000];
% percent_sig = [0.15, 0.20, 0.25];
percent_sig = [0.15, 0.2, 0.25];

for md = maxDemand_vec
    for psig = percent_sig
        gen_dem_data_50(md, psig);
    end
end
