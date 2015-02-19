maxDemand_vec = [1500, 2500, 3500, 4500];
percent_sig = [0.15, 0.20, 0.25];
% maxDemand_vec = [2000];
% percent_sig = [0.15];

for md = maxDemand_vec
    for psig = percent_sig
        gen_dem_data_15(md, psig);
    end
end
