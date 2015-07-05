maxDemand_vec = [3000,4000,5000,6000];
% maxDemand_vec = [3000];
percent_sig = [0.15, 0.20, 0.25];
% maxDemand_vec = [1500];
% percent_sig = [0.15];

for md = maxDemand_vec
    for psig = percent_sig
        gen_dem_data_30(md, psig);
    end
end
