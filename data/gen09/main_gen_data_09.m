maxDemand_vec = [2000, 3000, 4000, 5000];
percent_sig = [0.15, 0.20, 0.25];
% maxDemand_vec = [2000];
% percent_sig = [0.15];

for md = maxDemand_vec
    for psig = percent_sig
        gen_dem_data_09(md, psig);
    end
end
