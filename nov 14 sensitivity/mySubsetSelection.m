
p = 3;
combinations = combnk(1:37, p);
valid = zeros(size(combinations, 1),1);

for i=1:size(combinations, 1)
    % calculate rank
    % assume given dPas/dmu
    if rank(dPasdmu) == p
        valid(i) = 1;
    end
end