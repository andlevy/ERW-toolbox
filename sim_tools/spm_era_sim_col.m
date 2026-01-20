function col = spm_era_sim_col()


cold    = [
    142  69 163;
    141 182  20;
    242 102  20;
    20 149 170] / 255;
cold    = brighten(cold,-.1);
colb1   = brighten(cold,.66);
colb2   = brighten(colb1,.66);

col = [];
for c = 1:4

    col = [col;cold(c,:);colb1(c,:);colb2(c,:)];
end
